using System;
using System.Collections.Generic;
using System.Linq;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

namespace bertini_real
{
    /// <summary>
    /// Applies an ordered sequence of mesh boolean operations to each piece solid -- a fold,
    /// left to right.  Per piece: start from the Solid mesh, then for each Tool in order, union
    /// it (sign +1) or subtract it (sign <= 0) from the running result.  Order matters: each step
    /// acts on the result of the previous one (e.g. union body, then subtract the hole through it).
    ///
    /// Operations default to subtract, so feeding only negatives "just subtracts" them.  Booleans
    /// need closed solids: a non-closed Solid is a hard warning.  Features may be meshes or Breps
    /// (Breps are meshed first).  ("Feature" in the solid-modeling sense: an ordered additive or
    /// subtractive operation on a body -- a plug body is additive, a wire hole subtractive.)
    /// </summary>
    public class SurfaceBooleanPiece : GH_Component
    {
        public SurfaceBooleanPiece()
          : base("Boolean Piece", "BoolPiece",
                 "Fold an ordered list of mesh boolean ops (sign +1 union / -1 subtract) onto each piece solid",
                 "bertini_real", "Surface")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Solid", "S", "Closed piece mesh per piece (from Close Piece)", GH_ParamAccess.tree);
            pManager.AddGeometryParameter("Features", "F", "Geometry to boolean in, in order, per piece (meshes or Breps)", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("Operations", "O", "Sign per feature, parallel to Features: +1 union, -1 subtract. Defaults to subtract.", GH_ParamAccess.tree);
            Params.Input[1].Optional = true;   // no features -> pass the solid through unchanged
            Params.Input[2].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Result", "R", "Boolean result per piece", GH_ParamAccess.tree);
            pManager.AddBooleanParameter("Is Closed", "X", "Whether each piece's result is a closed solid", GH_ParamAccess.tree);
            pManager.AddTextParameter("Report", "Rep", "Per piece: any step that failed", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<GH_Mesh> solids;
            GH_Structure<IGH_GeometricGoo> features;
            GH_Structure<GH_Integer> ops;
            if (!DA.GetDataTree(0, out solids)) return;
            DA.GetDataTree(1, out features);  // optional; no features -> solid passes through
            DA.GetDataTree(2, out ops);       // optional; default sign is -1 (subtract)

            var resultTree = new DataTree<Mesh>();
            var closedTree = new DataTree<bool>();
            var reportTree = new DataTree<string>();

            for (int b = 0; b < solids.PathCount; b++)
            {
                GH_Path path = solids.get_Path(b);

                // starting solid(s)
                var current = new List<Mesh>();
                foreach (var goo in solids.get_Branch(path))
                    if (goo is GH_Mesh gm && gm.Value != null)
                        current.Add(gm.Value.DuplicateMesh());

                if (current.Count == 0) continue;
                if (current.Any(m => !m.IsClosed))
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                        $"Piece at path {path}: solid is not closed; boolean results are unreliable on open meshes.");

                // features + signs for this piece, in order
                var featureMeshes = (features != null && features.PathExists(path)) ? ToMeshList(features.get_Branch(path)) : new List<Mesh>();
                var signs = SignsForPath(ops, path, featureMeshes.Count);

                for (int i = 0; i < featureMeshes.Count; i++)
                {
                    Mesh feature = featureMeshes[i];
                    bool union = signs[i] > 0;

                    if (!feature.IsClosed)
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                            $"Piece at path {path}: feature {i} is not a closed mesh; the boolean may do nothing. " +
                            "Check the feature is a closed solid and actually overlaps the piece.");

                    Mesh[] next = union
                        ? Mesh.CreateBooleanUnion(current.Concat(new[] { feature }))
                        : Mesh.CreateBooleanDifference(current, new[] { feature });

                    if (next == null || next.Length == 0)
                    {
                        reportTree.Add($"step {i} ({(union ? "union" : "subtract")}) failed", path);
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                            $"Piece at path {path}: boolean step {i} ({(union ? "union" : "subtract")}) failed; kept the previous result.");
                        continue; // keep the last good `current`
                    }
                    current = next.ToList();
                }

                foreach (var m in current)
                    resultTree.Add(m, path);
                closedTree.Add(current.All(m => m != null && m.IsClosed), path);
            }

            DA.SetDataTree(0, resultTree);
            DA.SetDataTree(1, closedTree);
            DA.SetDataTree(2, reportTree);
        }

        /// <summary>Convert a branch of geometry goos to meshes (meshing any Breps).</summary>
        private static List<Mesh> ToMeshList(System.Collections.IList branch)
        {
            var meshes = new List<Mesh>();
            foreach (var goo in branch)
            {
                if (goo is GH_Mesh gm && gm.Value != null)
                {
                    meshes.Add(gm.Value);
                    continue;
                }
                GeometryBase geo = (goo as IGH_GeometricGoo)?.IsValid == true
                    ? GH_Convert.ToGeometryBase(goo)
                    : null;
                if (geo is Mesh m)
                {
                    meshes.Add(m);
                }
                else if (geo is Brep brep)
                {
                    Mesh[] fromBrep = Mesh.CreateFromBrep(brep, MeshingParameters.Default);
                    if (fromBrep != null && fromBrep.Length > 0)
                    {
                        // CreateFromBrep returns one mesh per face; appending leaves the cap
                        // seams unwelded (an open mesh), which makes mesh booleans no-op. Weld
                        // coincident vertices, and if it's still open (an uncapped tube, or a
                        // non-conforming seam between lateral + caps), fill the holes so the
                        // cutter becomes a closed solid.
                        var combined = new Mesh();
                        foreach (var mm in fromBrep) combined.Append(mm);
                        combined.Vertices.CombineIdentical(true, true);
                        if (!combined.IsClosed) combined.FillHoles();
                        combined.RebuildNormals();
                        combined.Compact();
                        meshes.Add(combined);
                    }
                }
            }
            return meshes;
        }

        /// <summary>Signs for a piece's tools; missing/short entries default to -1 (subtract).</summary>
        private static List<int> SignsForPath(GH_Structure<GH_Integer> ops, GH_Path path, int count)
        {
            var signs = new List<int>();
            System.Collections.IList branch = (ops != null && ops.PathExists(path)) ? ops.get_Branch(path) : null;
            for (int i = 0; i < count; i++)
            {
                int s = -1;
                if (branch != null && i < branch.Count && branch[i] is GH_Integer gi)
                    s = gi.Value;
                signs.Add(s);
            }
            return signs;
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("lego.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("6F2C1A94-8E3D-4B57-9A0C-2D7E5B41F8C3"); }
        }
    }
}
