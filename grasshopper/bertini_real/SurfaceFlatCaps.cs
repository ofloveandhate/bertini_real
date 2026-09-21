using System;
using System.Collections.Generic;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

namespace bertini_real
{
    /// <summary>
    /// Caps each surface piece's on-sphere boundary loop with a FLAT fan to the loop's centroid,
    /// rather than a cap that hugs the sphere (Sphere Caps).  The sphere itself still shows where
    /// the piece was cut; this just fills the opening flat -- so a Bertini-computed cylinder gets
    /// flat disk ends, and blocky decompositions get faceted flat caps.
    ///
    /// Same on-sphere boundary detection as Sphere Caps; the cap is built on the piece mesh's own
    /// boundary vertices (so it welds watertight in Close Piece), with the apex at the loop
    /// centroid.  Resolution adds concentric (linearly interpolated) rings for a denser flat cap.
    /// </summary>
    public class SurfaceFlatCaps : GH_Component
    {
        public SurfaceFlatCaps()
          : base("Flat Caps", "FlatCaps",
                 "Cap a surface piece's on-sphere boundary loops with a flat fan to the loop centroid",
                 "bertini_real", "Surface")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Meshes", "M", "Surface piece meshes (from Surface Read GH JSON)", GH_ParamAccess.tree);
            pManager.AddBrepParameter("Sphere", "S", "Bounding sphere Brep (from Surface Read GH JSON)", GH_ParamAccess.item);
            pManager.AddNumberParameter("Tolerance", "T", "Distance tolerance for testing whether a boundary vertex lies on the sphere", GH_ParamAccess.item, 1e-3);
            pManager.AddIntegerParameter("Resolution", "R", "Concentric rings from boundary to centroid (1 = a single flat fan)", GH_ParamAccess.item, 1);
            Params.Input[2].Optional = true;
            Params.Input[3].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Caps", "C", "Flat cap mesh(es) per piece (one per on-sphere loop)", GH_ParamAccess.tree);
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

            int resolution = 1;
            DA.GetData(3, ref resolution);
            if (resolution < 1) resolution = 1;

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

                foreach (var goo in meshes.get_Branch(path))
                {
                    var gm = goo as GH_Mesh;
                    if (gm?.Value == null) continue;
                    Mesh mesh = gm.Value;

                    List<List<int>> loops = Capping.OnSphereBoundaryLoops(mesh, center, radius, tol, out bool unclean);
                    if (unclean)
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                            $"Piece at path {path}: on-sphere boundary is not all clean closed loops; those parts were not capped.");

                    foreach (var loop in loops)
                    {
                        Mesh cap = BuildFlatCap(mesh, loop, resolution);
                        if (cap != null && cap.Faces.Count > 0)
                            caps.Add(cap, path);
                    }
                }
            }

            DA.SetDataTree(0, caps);
        }

        /// <summary>
        /// Flat fan over a closed loop: apex at the loop centroid, with <paramref name="resolution"/>
        /// concentric rings linearly interpolated from the (fixed, shared) boundary toward the
        /// centroid.  Ring 0 keeps the exact boundary vertices so a later join welds watertight.
        /// </summary>
        private static Mesh BuildFlatCap(Mesh mesh, List<int> loop, int resolution)
        {
            var topo = mesh.TopologyVertices;
            int n = loop.Count;

            var pts = new Point3d[n];
            double cx = 0, cy = 0, cz = 0;
            for (int k = 0; k < n; k++)
            {
                pts[k] = topo[loop[k]];
                cx += pts[k].X; cy += pts[k].Y; cz += pts[k].Z;
            }
            var centroid = new Point3d(cx / n, cy / n, cz / n);

            int R = Math.Max(1, resolution);
            var cap = new Mesh();

            // ring 0: exact boundary; rings 1..R-1: linearly interpolated toward the centroid
            for (int k = 0; k < n; k++) cap.Vertices.Add(pts[k]);
            for (int r = 1; r < R; r++)
            {
                double t = (double)r / R;
                for (int k = 0; k < n; k++)
                    cap.Vertices.Add(pts[k] + t * (centroid - pts[k]));
            }
            int apexIdx = cap.Vertices.Add(centroid);

            int Idx(int r, int k) => r * n + k;

            for (int r = 0; r < R - 1; r++)
            {
                for (int k = 0; k < n; k++)
                {
                    int k2 = (k + 1) % n;
                    Capping.AddTri(cap, Idx(r, k), Idx(r, k2), Idx(r + 1, k2));
                    Capping.AddTri(cap, Idx(r, k), Idx(r + 1, k2), Idx(r + 1, k));
                }
            }
            for (int k = 0; k < n; k++)
            {
                int k2 = (k + 1) % n;
                Capping.AddTri(cap, Idx(R - 1, k), Idx(R - 1, k2), apexIdx);
            }

            if (cap.Faces.Count == 0) return null;
            cap.Normals.ComputeNormals();
            cap.Compact();
            return cap;
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("lego.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("1D4A7E62-9C58-4B03-8E71-2F60A9C4D5B8"); }
        }
    }
}
