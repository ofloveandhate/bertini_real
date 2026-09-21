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
    /// decompose into constituents.)  Each piece is one tree branch and is translated radially
    /// outward from the common center by Distance model units (an absolute distance, matching
    /// Spread By Connectors): Distance = 0 leaves everything in place, larger spreads them further.
    ///
    /// Operates on any per-piece GEOMETRY tree, so it accepts bare meshes (Surface Read GH JSON /
    /// Close Piece) or a piece's mesh together with its connectors (Surface Group By Piece) -- all
    /// items in a branch move together, so a piece and its connectors stay assembled.
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
            pManager.AddGeometryParameter("Geometry", "G", "Per-piece geometry, one branch per piece: meshes, or a piece's mesh + connectors from Surface Group By Piece", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Distance", "D", "Distance each piece moves outward from the center (model units). 0 = no move.", GH_ParamAccess.item, 1.0);
            pManager.AddPointParameter("Center", "C", "Center to spread away from (default: average of the piece centers)", GH_ParamAccess.item);
            Params.Input[2].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGeometryParameter("Geometry", "G", "Spread-apart geometry (same tree structure)", GH_ParamAccess.tree);
            pManager.AddVectorParameter("Translations", "T", "Translation applied to each piece", GH_ParamAccess.tree);
            pManager.AddPointParameter("Center", "C", "The center the pieces were spread from", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<IGH_GeometricGoo> geometry;
            if (!DA.GetDataTree(0, out geometry)) return;

            double distance = 1.0;
            DA.GetData(1, ref distance);

            Point3d center = Point3d.Unset;
            bool hasCenter = DA.GetData(2, ref center);

            // per-piece center = center of the branch's combined bounding box (all geometry in it)
            var paths = new List<GH_Path>();
            var pieceCenter = new List<Point3d>();
            for (int b = 0; b < geometry.PathCount; b++)
            {
                GH_Path path = geometry.get_Path(b);
                BoundingBox bb = BoundingBox.Empty;
                bool any = false;
                foreach (var goo in geometry.get_Branch(path))
                {
                    if (!(goo is IGH_GeometricGoo gg) || !gg.IsValid) continue;
                    bb.Union(gg.Boundingbox);
                    any = true;
                }
                if (!any) continue;
                paths.Add(path);
                pieceCenter.Add(bb.Center);
            }

            if (paths.Count == 0) return;

            if (paths.Count == 1)
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                    "Only one tree branch received, so there is nothing to spread (one piece per branch). " +
                    "The per-piece tree structure was probably flattened upstream (Scale, a wire, or a Flatten on this input) -- " +
                    "keep one branch per piece coming in.");

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

            var outGeo = new DataTree<IGH_GeometricGoo>();
            var outVec = new DataTree<Vector3d>();

            for (int i = 0; i < paths.Count; i++)
            {
                GH_Path path = paths[i];
                // absolute distance: move each piece Distance units along its outward (unit) direction
                Vector3d dir = pieceCenter[i] - origin;
                if (!dir.IsTiny()) dir.Unitize();
                Vector3d t = distance * dir;
                Transform xf = Transform.Translation(t);

                foreach (var goo in geometry.get_Branch(path))
                {
                    if (!(goo is IGH_GeometricGoo gg) || !gg.IsValid) continue;
                    // duplicate so the input geometry is left untouched, then translate
                    IGH_GeometricGoo moved = gg.DuplicateGeometry();
                    moved = moved.Transform(xf);
                    outGeo.Add(moved, path);
                }

                outVec.Add(t, path);
            }

            DA.SetDataTree(0, outGeo);
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
