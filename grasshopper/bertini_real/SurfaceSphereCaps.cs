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

                    List<List<int>> loops = Capping.OnSphereBoundaryLoops(mesh, center, radius, tol, out bool unclean);
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

            // The loop is traversed in an arbitrary direction, so the natural winding below may
            // come out facing the sphere center (inward).  The cap closes a solid that sits INSIDE
            // the sphere, so its outward normal must point radially away from the center.  Reverse
            // the loop order when the winding would be inward, so the cap is outward by construction
            // and agrees with the (outward) piece -- no seam fold for the later UnifyNormals to fix.
            if (CapWindsInward(pts, pole, center))
            {
                Array.Reverse(pts);
                Array.Reverse(dirs);
            }

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
                    Capping.AddTri(cap, Idx(r, k), Idx(r, k2), Idx(r + 1, k2));
                    Capping.AddTri(cap, Idx(r, k), Idx(r + 1, k2), Idx(r + 1, k));
                }
            }
            // innermost ring fans to the pole
            for (int k = 0; k < n; k++)
            {
                int k2 = (k + 1) % n;
                Capping.AddTri(cap, Idx(R - 1, k), Idx(R - 1, k2), poleIdx);
            }

            if (cap.Faces.Count == 0) return null;
            cap.Normals.ComputeNormals();
            cap.Compact();
            return cap;
        }

        /// <summary>
        /// True when the cap built from <paramref name="pts"/> (in their current order, fanned toward
        /// <paramref name="pole"/>) would have its faces pointing toward the sphere center instead of
        /// away from it.  The actual cap triangles share the orientation of the simple cone fan
        /// (pts[k] -> pts[k+1] -> pole), so we sum that fan's face normals dotted with the outward
        /// radial direction; a negative total means the winding is inward and the loop should reverse.
        /// </summary>
        private static bool CapWindsInward(Point3d[] pts, Point3d pole, Point3d center)
        {
            int n = pts.Length;
            double radialDot = 0.0;
            for (int k = 0; k < n; k++)
            {
                int k2 = (k + 1) % n;
                Vector3d nrm = Vector3d.CrossProduct(pts[k2] - pts[k], pole - pts[k]);
                Point3d mid = 0.5 * (pts[k] + pts[k2]);   // edge midpoint
                Vector3d outward = mid - center;          // radially outward from the sphere center
                radialDot += nrm * outward;
            }
            return radialDot < 0.0;
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
