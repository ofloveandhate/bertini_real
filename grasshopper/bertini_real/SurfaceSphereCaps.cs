using System;
using System.Collections.Generic;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino;
using Rhino.Geometry;

namespace bertini_real
{
    /// <summary>
    /// Builds the spherical "cap" mesh(es) that close a surface piece where it meets the
    /// decomposition's bounding sphere.  The cap is built on the piece mesh's OWN naked boundary
    /// loop that lies on the sphere -- never on the separately-sampled sphere curve -- so a raw
    /// piece gets a raw cap and a sampled piece a sampled cap (raw with raw, sampled with sampled),
    /// and the boundary vertices are shared so a later join welds watertight.
    ///
    /// For each closed on-sphere boundary loop, two candidate fan caps (one to each pole of the
    /// loop's mean direction) are built and the smaller-area one is kept.  Works in mesh topology
    /// space so coincident vertices (nodal singularities) merge, and skips degenerate edges/faces.
    /// </summary>
    public class SurfaceSphereCaps : GH_Component
    {
        public SurfaceSphereCaps()
          : base("Sphere Caps", "SphereCaps",
                 "Cap a surface piece's on-sphere boundary loops with faceted spherical caps (smaller-area side)",
                 "bertini_real", "Surface")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Meshes", "M", "Surface piece meshes (from Surface Read GH JSON)", GH_ParamAccess.tree);
            pManager.AddBrepParameter("Sphere", "S", "Bounding sphere Brep (from Surface Read GH JSON)", GH_ParamAccess.item);
            pManager.AddNumberParameter("Tolerance", "T", "Distance tolerance for testing whether a boundary vertex lies on the sphere", GH_ParamAccess.item, 1e-3);
            Params.Input[2].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Caps", "C", "Spherical cap mesh(es) per piece (one per on-sphere loop)", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<GH_Mesh> meshes;
            if (!DA.GetDataTree(0, out meshes)) return;

            Brep sphereBrep = null;
            if (!DA.GetData(1, ref sphereBrep) || sphereBrep == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "No sphere supplied.");
                return;
            }

            double tol = 1e-3;
            DA.GetData(2, ref tol);

            // recover center/radius from the sphere Brep's bounding box (exact for a sphere)
            BoundingBox bb = sphereBrep.GetBoundingBox(true);
            Point3d center = bb.Center;
            double radius = bb.Diagonal.X / 2.0;
            if (radius <= 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Sphere has non-positive radius.");
                return;
            }

            var caps = new DataTree<Mesh>();

            for (int b = 0; b < meshes.PathCount; b++)
            {
                GH_Path path = meshes.get_Path(b);
                var branch = meshes.get_Branch(path);

                foreach (var goo in branch)
                {
                    var gm = goo as GH_Mesh;
                    if (gm?.Value == null) continue;
                    Mesh mesh = gm.Value;

                    List<List<int>> loops = OnSphereBoundaryLoops(mesh, center, radius, tol, out bool unclean);
                    if (unclean)
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                            $"Piece at path {path}: on-sphere boundary is not all clean closed loops (open arc or non-manifold junction); those parts were not capped.");

                    foreach (var loop in loops)
                    {
                        Mesh cap = BuildSmallerCap(mesh, loop, center, radius, tol);
                        if (cap != null && cap.Faces.Count > 0)
                            caps.Add(cap, path);
                    }
                }
            }

            DA.SetDataTree(0, caps);
        }

        /// <summary>
        /// Returns the closed loops (as ordered lists of TopologyVertex indices) of the mesh's
        /// naked boundary that lie on the sphere.  Sets <paramref name="unclean"/> if some on-sphere
        /// boundary did not form clean degree-2 cycles.  Works in topology space (coincident
        /// vertices merged); zero-length edges collapse and are ignored.
        /// </summary>
        private static List<List<int>> OnSphereBoundaryLoops(Mesh mesh, Point3d center, double radius, double tol, out bool unclean)
        {
            unclean = false;
            var result = new List<List<int>>();

            var topo = mesh.TopologyVertices;
            var edges = mesh.TopologyEdges;

            bool OnSphere(int tv)
            {
                Point3d p = topo[tv];
                return Math.Abs(p.DistanceTo(center) - radius) < tol;
            }

            // adjacency among on-sphere naked topology edges
            var adj = new Dictionary<int, List<int>>();
            void Link(int u, int v)
            {
                if (!adj.TryGetValue(u, out var lu)) { lu = new List<int>(); adj[u] = lu; }
                if (!lu.Contains(v)) lu.Add(v);
            }

            for (int e = 0; e < edges.Count; e++)
            {
                if (edges.GetConnectedFaces(e).Length != 1) continue; // not naked
                IndexPair ip = edges.GetTopologyVertices(e);
                int i = ip.I, j = ip.J;
                if (i == j) continue;                 // degenerate (collapsed) edge
                if (!OnSphere(i) || !OnSphere(j)) continue;
                Link(i, j);
                Link(j, i);
            }

            if (adj.Count == 0) return result;

            // walk connected components into ordered cycles
            var visited = new HashSet<int>();
            foreach (int start in adj.Keys)
            {
                if (visited.Contains(start)) continue;

                var loop = new List<int>();
                int prev = -1, cur = start;
                bool clean = true;

                while (true)
                {
                    visited.Add(cur);
                    loop.Add(cur);

                    var nbrs = adj[cur];
                    if (nbrs.Count != 2) { clean = false; break; } // junction or dead-end

                    int next = nbrs[0] != prev ? nbrs[0] : nbrs[1];
                    if (next == start) break;          // closed the loop
                    if (visited.Contains(next)) { clean = false; break; }
                    prev = cur;
                    cur = next;
                }

                if (clean && loop.Count >= 3)
                    result.Add(loop);
                else
                    unclean = true;
            }

            return result;
        }

        /// <summary>
        /// Builds the smaller-area fan cap over a closed loop of topology-vertex indices, with the
        /// apex on the sphere.  Skips degenerate (near-zero-area) triangles.
        /// </summary>
        private static Mesh BuildSmallerCap(Mesh mesh, List<int> loop, Point3d center, double radius, double tol)
        {
            var topo = mesh.TopologyVertices;
            int n = loop.Count;

            var pts = new Point3d[n];
            var mean = new Vector3d(0, 0, 0);
            for (int k = 0; k < n; k++)
            {
                pts[k] = topo[loop[k]];
                mean += pts[k] - center;
            }
            if (mean.IsTiny()) mean = new Vector3d(0, 0, 1);
            mean.Unitize();

            Mesh best = null;
            double bestArea = double.MaxValue;

            foreach (int sign in new[] { 1, -1 })
            {
                Point3d apex = center + sign * radius * mean;

                var cap = new Mesh();
                cap.Vertices.AddVertices(pts);
                int apexIdx = cap.Vertices.Add(apex);

                double area = 0.0;
                for (int k = 0; k < n; k++)
                {
                    int a = k, c = (k + 1) % n;
                    double triArea = 0.5 * Vector3d.CrossProduct(pts[a] - apex, pts[c] - apex).Length;
                    if (triArea < tol * tol) continue;   // degenerate triangle -> skip
                    cap.Faces.AddFace(a, c, apexIdx);
                    area += triArea;
                }

                if (cap.Faces.Count == 0) continue;

                if (area < bestArea)
                {
                    bestArea = area;
                    best = cap;
                }
            }

            if (best != null)
            {
                best.Normals.ComputeNormals();
                best.Compact();
            }
            return best;
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("lego.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("A7D2E4F1-6C90-4B33-9E27-5C8F1A0B7E64"); }
        }
    }
}
