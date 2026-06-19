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
    /// Spreads the pieces of a decomposition apart so they can be seen individually without
    /// touching.  (Named "Spread" rather than "Explode" to avoid Grasshopper's sense of explode =
    /// decompose into constituents.)  Each piece (one tree branch) is translated radially away
    /// from the common center by Factor times its offset from that center: Factor = 0 leaves
    /// everything in place, larger Factor spreads them further.  All meshes in a branch move
    /// together, so a piece stays intact.
    /// </summary>
    public class SurfaceSpreadPieces : GH_Component
    {
        public SurfaceSpreadPieces()
          : base("Spread Pieces", "Spread",
                 "Move decomposition pieces apart (radially from their common center) so they can be seen separated",
                 "bertini_real", "Surface")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Meshes", "M", "Per-piece meshes (one branch per piece)", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Factor", "F", "Spread amount: each piece moves by Factor x (its center - the overall center). 0 = no move.", GH_ParamAccess.item, 0.5);
            pManager.AddPointParameter("Center", "C", "Center to spread away from (default: average of the piece centers)", GH_ParamAccess.item);
            Params.Input[2].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Meshes", "M", "Spread-apart piece meshes (same tree structure)", GH_ParamAccess.tree);
            pManager.AddVectorParameter("Translations", "T", "Translation applied to each piece", GH_ParamAccess.tree);
            pManager.AddPointParameter("Center", "C", "The center the pieces were spread from", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<GH_Mesh> meshes;
            if (!DA.GetDataTree(0, out meshes)) return;

            double factor = 0.5;
            DA.GetData(1, ref factor);

            Point3d center = Point3d.Unset;
            bool hasCenter = DA.GetData(2, ref center);

            // per-branch (per-piece) center = center of the branch's combined bounding box
            var paths = new List<GH_Path>();
            var pieceCenter = new List<Point3d>();
            for (int b = 0; b < meshes.PathCount; b++)
            {
                GH_Path path = meshes.get_Path(b);
                BoundingBox bb = BoundingBox.Empty;
                bool any = false;
                foreach (var goo in meshes.get_Branch(path))
                    if (goo is GH_Mesh gm && gm.Value != null)
                    {
                        bb.Union(gm.Value.GetBoundingBox(false));
                        any = true;
                    }
                if (!any) continue;
                paths.Add(path);
                pieceCenter.Add(bb.Center);
            }

            if (paths.Count == 0) return;

            // overall center: explicit input, else the average of the piece centers
            Point3d origin;
            if (hasCenter)
            {
                origin = center;
            }
            else
            {
                double x = 0, y = 0, z = 0;
                foreach (var c in pieceCenter) { x += c.X; y += c.Y; z += c.Z; }
                origin = new Point3d(x / pieceCenter.Count, y / pieceCenter.Count, z / pieceCenter.Count);
            }

            var outMesh = new DataTree<Mesh>();
            var outVec = new DataTree<Vector3d>();

            for (int i = 0; i < paths.Count; i++)
            {
                GH_Path path = paths[i];
                Vector3d t = factor * (pieceCenter[i] - origin);
                Transform xf = Transform.Translation(t);

                foreach (var goo in meshes.get_Branch(path))
                    if (goo is GH_Mesh gm && gm.Value != null)
                    {
                        Mesh m = gm.Value.DuplicateMesh();
                        m.Transform(xf);
                        outMesh.Add(m, path);
                    }

                outVec.Add(t, path);
            }

            DA.SetDataTree(0, outMesh);
            DA.SetDataTree(1, outVec);
            DA.SetData(2, origin);
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("transform.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("C2A6F39B-4E78-4D15-8B0A-3F9E1C7D6052"); }
        }
    }
}
