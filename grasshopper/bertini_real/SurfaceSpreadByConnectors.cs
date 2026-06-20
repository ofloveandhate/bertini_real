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
    /// Spreads pieces along their connector axes, to show how they would assemble -- an exploded
    /// *assembly* view, as opposed to the purely radial Spread Pieces.
    ///
    /// For each piece, every incident singularity contributes a unit vector along its connector
    /// direction, oriented AWAY from the singularity (the way the piece slides off its rod); the
    /// piece is translated by Factor times the sum.  A piece connected on opposite sides (e.g. a
    /// central hub) barely moves because its axes cancel, while leaf pieces slide out along their
    /// rods.  Wire Sing Locations / Sing Directions / Sing On Pieces straight from the reader.
    /// </summary>
    public class SurfaceSpreadByConnectors : GH_Component
    {
        public SurfaceSpreadByConnectors()
          : base("Spread By Connectors", "SpreadConn",
                 "Spread pieces along their connector axes (assembly explosion) using the singularity directions",
                 "bertini_real", "Surface")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGeometryParameter("Geometry", "G", "Per-piece geometry, one branch per piece (meshes, or mesh + connectors)", GH_ParamAccess.tree);
            pManager.AddPointParameter("Sing Locations", "SL", "Singularity locations (from Surface Read GH JSON)", GH_ParamAccess.list);
            pManager.AddVectorParameter("Sing Directions", "SD", "Singularity connector directions (from Surface Read GH JSON)", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Sing On Pieces", "SOP", "Per piece: indices of the singularities on it (from Surface Read GH JSON)", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Factor", "F", "Separation distance per connector (model units). 0 = no move.", GH_ParamAccess.item, 1.0);
            Params.Input[4].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGeometryParameter("Geometry", "G", "Spread-apart geometry (same tree structure)", GH_ParamAccess.tree);
            pManager.AddVectorParameter("Translations", "T", "Translation applied to each piece", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<IGH_GeometricGoo> geometry;
            if (!DA.GetDataTree(0, out geometry)) return;

            var locations = new List<Point3d>();
            var directions = new List<Vector3d>();
            GH_Structure<GH_Integer> onPieces;
            DA.GetDataList(1, locations);
            DA.GetDataList(2, directions);
            if (!DA.GetDataTree(3, out onPieces)) return;

            double factor = 1.0;
            DA.GetData(4, ref factor);

            var outGeo = new DataTree<IGH_GeometricGoo>();
            var outVec = new DataTree<Vector3d>();

            for (int b = 0; b < geometry.PathCount; b++)
            {
                GH_Path path = geometry.get_Path(b);

                // piece center (combined bounding box of its geometry)
                BoundingBox bb = BoundingBox.Empty;
                bool any = false;
                foreach (var goo in geometry.get_Branch(path))
                    if (goo is IGH_GeometricGoo gg && gg.IsValid) { bb.Union(gg.Boundingbox); any = true; }
                if (!any) continue;
                Point3d centroid = bb.Center;

                // sum a unit vector per incident singularity, oriented away from the singularity
                Vector3d disp = Vector3d.Zero;
                if (onPieces.PathExists(path))
                {
                    foreach (var goo in onPieces.get_Branch(path))
                    {
                        if (!(goo is GH_Integer gi)) continue;
                        int idx = gi.Value;
                        if (idx < 0 || idx >= locations.Count || idx >= directions.Count) continue;

                        Vector3d axis = directions[idx];
                        if (axis.IsTiny()) continue;
                        axis.Unitize();

                        // orient the axis so the piece moves away from the singularity
                        double along = (centroid - locations[idx]) * axis;
                        disp += (along >= 0 ? 1.0 : -1.0) * axis;
                    }
                }

                disp *= factor;
                Transform xf = Transform.Translation(disp);

                foreach (var goo in geometry.get_Branch(path))
                {
                    if (!(goo is IGH_GeometricGoo gg) || !gg.IsValid) continue;
                    IGH_GeometricGoo moved = gg.DuplicateGeometry();
                    moved = moved.Transform(xf);
                    outGeo.Add(moved, path);
                }
                outVec.Add(disp, path);
            }

            DA.SetDataTree(0, outGeo);
            DA.SetDataTree(1, outVec);
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("transform.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("9F26C7B4-3A81-4D60-8E15-7C0B2F95E6A3"); }
        }
    }
}
