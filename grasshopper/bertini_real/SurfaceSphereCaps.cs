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
            pManager.AddIntegerParameter("Resolution", "R", "Radial subdivisions of the cap (rings from boundary to pole); higher = smoother, follows the sphere", GH_ParamAccess.item, 4);
            Params.Input[2].Optional = true;
            Params.Input[3].Optional = true;
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

            int resolution = 4;
            DA.GetData(3, ref resolution);
            if (resolution < 1) resolution = 1;

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
                        Mesh cap = BuildSmallerCap(mesh, loop, center, radius, tol, resolution);
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
        /// Builds the smaller-area cap over a closed loop of topology-vertex indices.  The pole is
        /// placed on the sphere along the loop's mean direction; the cap is subdivided into
        /// <paramref name="resolution"/> radial rings, each slerped along the sphere from the
        /// (fixed, shared) boundary toward the pole, so it follows the sphere instead of pinching
        /// to a flat point.  Ring 0 keeps the exact boundary vertices so a later join welds
        /// watertight.  Degenerate (near-zero-area) triangles are skipped.
        /// </summary>
        private static Mesh BuildSmallerCap(Mesh mesh, List<int> loop, Point3d center, double radius, double tol, int resolution)
        {
            var topo = mesh.TopologyVertices;
            int n = loop.Count;

            var pts = new Point3d[n];      // exact boundary points (ring 0)
            var dirs = new Vector3d[n];    // unit directions from center
            var mean = new Vector3d(0, 0, 0);
            for (int k = 0; k < n; k++)
            {
                pts[k] = topo[loop[k]];
                Vector3d d = pts[k] - center;
                dirs[k] = d;
                dirs[k].Unitize();
                mean += d;
            }
            if (mean.IsTiny()) mean = new Vector3d(0, 0, 1);
            mean.Unitize();

            // pick the pole side that yields the smaller cap
            int sign = SmallerCapSign(pts, center, radius, mean);
            Vector3d poleDir = sign * mean;
            Point3d pole = center + radius * poleDir;

            int R = Math.Max(1, resolution);
            var cap = new Mesh();

            // ring 0: exact boundary; rings 1..R-1: slerped toward the pole
            for (int k = 0; k < n; k++) cap.Vertices.Add(pts[k]);
            for (int r = 1; r < R; r++)
            {
                double t = (double)r / R;
                for (int k = 0; k < n; k++)
                {
                    Vector3d dir = Slerp(dirs[k], poleDir, t);
                    cap.Vertices.Add(center + radius * dir);
                }
            }
            int poleIdx = cap.Vertices.Add(pole);

            int Idx(int r, int k) => r * n + k;

            // quad strips between consecutive full rings
            for (int r = 0; r < R - 1; r++)
            {
                for (int k = 0; k < n; k++)
                {
                    int k2 = (k + 1) % n;
                    AddTri(cap, Idx(r, k), Idx(r, k2), Idx(r + 1, k2));
                    AddTri(cap, Idx(r, k), Idx(r + 1, k2), Idx(r + 1, k));
                }
            }
            // innermost ring fans to the pole
            for (int k = 0; k < n; k++)
            {
                int k2 = (k + 1) % n;
                AddTri(cap, Idx(R - 1, k), Idx(R - 1, k2), poleIdx);
            }

            if (cap.Faces.Count == 0) return null;
            cap.Normals.ComputeNormals();
            cap.Compact();
            return cap;
        }

        private static int SmallerCapSign(Point3d[] pts, Point3d center, double radius, Vector3d mean)
        {
            double Area(int sign)
            {
                Point3d apex = center + sign * radius * mean;
                double area = 0.0;
                int n = pts.Length;
                for (int k = 0; k < n; k++)
                {
                    int k2 = (k + 1) % n;
                    area += 0.5 * Vector3d.CrossProduct(pts[k] - apex, pts[k2] - apex).Length;
                }
                return area;
            }
            return Area(1) <= Area(-1) ? 1 : -1;
        }

        // a tiny absolute distance below which two cap vertices are treated as the same point
        private const double CoincidentTol = 1e-9;

        private static void AddTri(Mesh m, int i, int j, int k)
        {
            Point3d a = m.Vertices[i];
            Point3d b = m.Vertices[j];
            Point3d c = m.Vertices[k];
            // skip only TRULY degenerate triangles (a collapsed edge); thin-but-valid triangles
            // must be kept, or high-resolution caps develop tiny holes near the boundary.
            if (a.DistanceTo(b) < CoincidentTol ||
                b.DistanceTo(c) < CoincidentTol ||
                a.DistanceTo(c) < CoincidentTol)
                return;
            m.Faces.AddFace(i, j, k);
        }

        /// <summary>Spherical interpolation of two unit vectors; linear fallback when (anti)parallel.</summary>
        private static Vector3d Slerp(Vector3d v0, Vector3d v1, double t)
        {
            double dot = Math.Max(-1.0, Math.Min(1.0, v0 * v1));
            double omega = Math.Acos(dot);
            double so = Math.Sin(omega);
            if (omega < 1e-9 || so < 1e-9)
            {
                Vector3d lin = (1.0 - t) * v0 + t * v1;
                if (lin.IsTiny()) return v0;
                lin.Unitize();
                return lin;
            }
            return (Math.Sin((1.0 - t) * omega) / so) * v0 + (Math.Sin(t * omega) / so) * v1;
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("lego.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("A7D2E4F1-6C90-4B33-9E27-5C8F1A0B7E64"); }
        }
    }
}
