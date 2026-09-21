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
    /// Joins a surface piece mesh with its spherical cap mesh(es) into a single welded mesh,
    /// merging coincident boundary vertices so the result is watertight where the cap meets the
    /// piece.  Reports whether each joined piece is a closed solid.  (A piece bounded only by the
    /// sphere closes fully; one that also abuts a singular curve stays open there until joined to
    /// its neighbor -- so IsClosed tells you which pieces are complete solids.)
    ///
    /// The piece meshes and caps must be the same sampling (raw with raw, sampled with sampled);
    /// since the caps are built from the piece's own boundary, feeding matching trees preserves that.
    /// </summary>
    public class SurfaceClosePiece : GH_Component
    {
        public SurfaceClosePiece()
          : base("Close Piece", "ClosePiece",
                 "Join a surface piece with its sphere cap(s) into a welded, ideally closed, mesh",
                 "bertini_real", "Surface")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Meshes", "M", "Surface piece meshes (from Surface Read GH JSON)", GH_ParamAccess.tree);
            pManager.AddMeshParameter("Caps", "C", "Sphere cap meshes per piece (from Sphere Caps)", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Tolerance", "T", "Vertex merge tolerance for welding cap to piece", GH_ParamAccess.item, 1e-6);
            Params.Input[2].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Closed", "M", "Piece joined with its cap(s), welded", GH_ParamAccess.tree);
            pManager.AddBooleanParameter("Is Closed", "X", "Whether the joined mesh is a closed solid", GH_ParamAccess.tree);
            pManager.AddCurveParameter("Naked Edges", "N", "Remaining naked (unwelded/open) edges, for diagnosing why a piece isn't closed", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<GH_Mesh> pieces;
            GH_Structure<GH_Mesh> caps;
            if (!DA.GetDataTree(0, out pieces)) return;
            DA.GetDataTree(1, out caps);   // caps may legitimately be empty for some pieces

            double tol = 1e-6;
            DA.GetData(2, ref tol);

            var outMesh = new DataTree<Mesh>();
            var outClosed = new DataTree<bool>();
            var outNaked = new DataTree<Polyline>();

            for (int b = 0; b < pieces.PathCount; b++)
            {
                GH_Path path = pieces.get_Path(b);

                var combined = new Mesh();

                foreach (var goo in pieces.get_Branch(path))
                    if (goo is GH_Mesh gm && gm.Value != null)
                        combined.Append(gm.Value);

                if (caps.PathExists(path))
                    foreach (var goo in caps.get_Branch(path))
                        if (goo is GH_Mesh gc && gc.Value != null)
                            combined.Append(gc.Value);

                if (combined.Faces.Count == 0)
                    continue;

                // weld: merge coincident vertices (the shared cap/piece boundary), drop junk
                combined.Vertices.CombineIdentical(true, true);
                combined.Faces.CullDegenerateFaces();
                combined.Vertices.CullUnused();
                combined.Compact();
                combined.RebuildNormals();
                combined.UnifyNormals();

                // UnifyNormals makes the faces mutually consistent but seeds from an arbitrary face,
                // so a closed solid can come out uniformly inside-out.  Mesh.Volume() is signed by
                // normal orientation, so a negative volume means the normals point inward -- flip the
                // whole mesh outward (matches the Python pipeline's trimesh.fix_normals()).
                bool closed = combined.IsClosed;
                if (closed && combined.Volume() < 0.0)
                {
                    combined.Flip(true, true, true);
                    combined.RebuildNormals();
                }
                Polyline[] naked = combined.GetNakedEdges() ?? Array.Empty<Polyline>();
                if (!closed)
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Remark,
                        $"Piece at path {path} is not a closed solid: {naked.Length} naked edge loop(s) remain " +
                        "(see Naked Edges output; likely an unwelded cap seam, an uncapped singular boundary, or a missing cap).");

                outMesh.Add(combined, path);
                outClosed.Add(closed, path);
                foreach (var pl in naked)
                    outNaked.Add(pl, path);
            }

            DA.SetDataTree(0, outMesh);
            DA.SetDataTree(1, outClosed);
            DA.SetDataTree(2, outNaked);
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("lego.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("D9C4B6A8-2E51-4F70-A38C-6B1D90E5F273"); }
        }
    }
}
