using System;
using System.Collections.Generic;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

namespace bertini_real
{
    /// <summary>
    /// Weaves the four connector trees from Surface Place Components into a single ordered
    /// Features tree plus matching Operations signs, ready to drop into Boolean Piece -- so you
    /// don't have to Merge/Weave/sign things by hand.
    ///
    /// Order, per piece: each connector as (positive +1, then negative -1), plugs then sockets --
    /// the pos, neg, pos, neg sequence.  Wire only the negatives and you get an all-subtract
    /// Features list (the short-term case).  Output Features/Operations are parallel and
    /// per-piece ({piece}).
    /// </summary>
    public class SurfaceConnectorsToFeatures : GH_Component
    {
        public SurfaceConnectorsToFeatures()
          : base("Connectors To Features", "Conn2Feat",
                 "Weave plug/socket connectors into one ordered Features tree + Operations signs for Boolean Piece",
                 "bertini_real", "Surface")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGeometryParameter("Plugs positive", "Plugs+", "Positive plugs per piece (union)", GH_ParamAccess.tree);
            pManager.AddGeometryParameter("Plugs negative", "Plugs-", "Negative plugs per piece (subtract)", GH_ParamAccess.tree);
            pManager.AddGeometryParameter("Sockets positive", "Sockets+", "Positive sockets per piece (union)", GH_ParamAccess.tree);
            pManager.AddGeometryParameter("Sockets negative", "Sockets-", "Negative sockets per piece (subtract)", GH_ParamAccess.tree);
            Params.Input[0].Optional = true;
            Params.Input[1].Optional = true;
            Params.Input[2].Optional = true;
            Params.Input[3].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGeometryParameter("Features", "F", "Ordered boolean features per piece (for Boolean Piece)", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("Operations", "O", "Sign per feature, parallel to Features: +1 union, -1 subtract", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<IGH_GeometricGoo> plugPos, plugNeg, socketPos, socketNeg;
            DA.GetDataTree(0, out plugPos);
            DA.GetDataTree(1, out plugNeg);
            DA.GetDataTree(2, out socketPos);
            DA.GetDataTree(3, out socketNeg);

            // union of piece indices present across all four inputs
            var pieceIndices = new SortedSet<int>();
            foreach (var tree in new[] { plugPos, plugNeg, socketPos, socketNeg })
                CollectPieceIndices(tree, pieceIndices);

            var featuresTree = new DataTree<IGH_GeometricGoo>();
            var opsTree = new DataTree<int>();

            foreach (int pieceIndex in pieceIndices)
            {
                GH_Path path = new GH_Path(pieceIndex);

                var pp = BranchItems(plugPos, path);
                var pn = BranchItems(plugNeg, path);
                var sp = BranchItems(socketPos, path);
                var sn = BranchItems(socketNeg, path);

                // plugs: (positive +1, negative -1) per connector, then sockets the same way
                Weave(pp, pn, path, featuresTree, opsTree);
                Weave(sp, sn, path, featuresTree, opsTree);
            }

            DA.SetDataTree(0, featuresTree);
            DA.SetDataTree(1, opsTree);
        }

        /// <summary>Emit (positive +1, negative -1) pairs in order; extras (unequal counts) follow with their own sign.</summary>
        private static void Weave(List<IGH_GeometricGoo> positives, List<IGH_GeometricGoo> negatives,
                                  GH_Path path, DataTree<IGH_GeometricGoo> features, DataTree<int> ops)
        {
            int n = Math.Max(positives.Count, negatives.Count);
            for (int i = 0; i < n; i++)
            {
                if (i < positives.Count) { features.Add(positives[i], path); ops.Add(+1, path); }
                if (i < negatives.Count) { features.Add(negatives[i], path); ops.Add(-1, path); }
            }
        }

        private static List<IGH_GeometricGoo> BranchItems(GH_Structure<IGH_GeometricGoo> tree, GH_Path path)
        {
            var items = new List<IGH_GeometricGoo>();
            if (tree != null && tree.PathExists(path))
                foreach (var goo in tree.get_Branch(path))
                    if (goo is IGH_GeometricGoo gg && gg.IsValid)
                        items.Add(gg);
            return items;
        }

        private static void CollectPieceIndices(GH_Structure<IGH_GeometricGoo> tree, SortedSet<int> into)
        {
            if (tree == null) return;
            for (int b = 0; b < tree.PathCount; b++)
            {
                var idx = tree.get_Path(b).Indices;
                into.Add(idx.Length > 0 ? idx[idx.Length - 1] : b);
            }
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("transform.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("4B8E2D17-5C6A-49F3-A1B0-3E9D7C625A48"); }
        }
    }
}
