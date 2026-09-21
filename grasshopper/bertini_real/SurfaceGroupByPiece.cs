using System;
using System.Collections;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

namespace bertini_real
{
    /// <summary>
    /// Groups each piece's mesh with its connectors into one branch per piece.  Everything
    /// upstream is already keyed per piece by tree branch (Surface Read GH JSON meshes,
    /// Surface Place Components connectors), so this just merges those trees branch-by-branch
    /// -- no JSON file, no pieceID-string routing.  Each output branch is: the mesh(es) for that
    /// piece followed by its connectors.
    /// </summary>
    public class SurfaceGroupByPiece : GH_Component
    {
        public SurfaceGroupByPiece()
          : base("Surface Group By Piece", "SurfGroupPiece",
              "Group each piece's mesh with its connectors (one branch per piece)",
              "bertini_real", "Surface")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGeometryParameter("Meshes", "M", "Piece meshes, one branch per piece (from Surface Read GH JSON)", GH_ParamAccess.tree);
            pManager.AddGeometryParameter("Plugs positive", "Plugs+", "Positive plugs per piece (from Surface Place Components)", GH_ParamAccess.tree);
            pManager.AddGeometryParameter("Plugs negative", "Plugs-", "Negative plugs per piece", GH_ParamAccess.tree);
            pManager.AddGeometryParameter("Sockets positive", "Sockets+", "Positive sockets per piece", GH_ParamAccess.tree);
            pManager.AddGeometryParameter("Sockets negative", "Sockets-", "Negative sockets per piece", GH_ParamAccess.tree);

            Params.Input[1].Optional = true;
            Params.Input[2].Optional = true;
            Params.Input[3].Optional = true;
            Params.Input[4].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGeometryParameter("Pieces with connectors", "PCs",
                "DataTree, one branch per piece: the mesh(es) followed by that piece's connectors.",
                GH_ParamAccess.tree);
            pManager.AddIntegerParameter("Piece Indices", "Is", "Piece index of each output branch, in branch order", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<IGH_GeometricGoo> meshes, plugsPos, plugsNeg, socketsPos, socketsNeg;
            if (!DA.GetDataTree(0, out meshes)) return;
            DA.GetDataTree(1, out plugsPos);
            DA.GetDataTree(2, out plugsNeg);
            DA.GetDataTree(3, out socketsPos);
            DA.GetDataTree(4, out socketsNeg);

            // piece index -> geometry (mesh first because meshes are ingested first)
            var grouped = new SortedDictionary<int, List<IGH_GeometricGoo>>();

            void ingest(GH_Structure<IGH_GeometricGoo> tree)
            {
                if (tree == null) return;
                for (int b = 0; b < tree.PathCount; b++)
                {
                    GH_Path path = tree.get_Path(b);
                    int pieceIndex = path.Indices.Length > 0 ? path.Indices[path.Indices.Length - 1] : b;

                    if (!grouped.TryGetValue(pieceIndex, out var items))
                    {
                        items = new List<IGH_GeometricGoo>();
                        grouped[pieceIndex] = items;
                    }
                    foreach (var goo in tree.get_Branch(path))
                        if (goo is IGH_GeometricGoo gg)
                            items.Add(gg);
                }
            }

            ingest(meshes);
            ingest(plugsPos);
            ingest(plugsNeg);
            ingest(socketsPos);
            ingest(socketsNeg);

            var outTree = new GH_Structure<IGH_GeometricGoo>();
            var pieceIndices = new List<int>();
            foreach (var kv in grouped)
            {
                GH_Path path = new GH_Path(kv.Key);
                foreach (var goo in kv.Value)
                    outTree.Append(goo, path);
                pieceIndices.Add(kv.Key);
            }

            DA.SetDataTree(0, outTree);
            DA.SetDataList(1, pieceIndices);
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("lego.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("9094345E-E59C-4E5D-B260-E24CFC11EFC3"); }
        }
    }
}
