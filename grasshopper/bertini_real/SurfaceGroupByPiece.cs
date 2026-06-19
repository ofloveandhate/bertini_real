using System;
using System.Collections.Generic;
using System.IO;
using System.Text.Json;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

namespace bertini_real
{
    public class SurfaceGroupByPiece : GH_Component
    {
        public SurfaceGroupByPiece()
          : base("Surface Group By Piece", "SurfGroupPiece",
              "Groups piece meshes with their transformed connectors",
              "bertini_real", "Surface")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("File Path", "F", "Path of json file with surface specs", GH_ParamAccess.item);
            pManager.AddGeometryParameter("Meshes", "Ms", "Piece meshes in piece order from SurfaceImportStlPieces", GH_ParamAccess.list);
            pManager.AddGeometryParameter("Plugs positive", "Plugs+", "Transformed positive plugs", GH_ParamAccess.list);
            pManager.AddGeometryParameter("Plugs negative", "Plugs-", "Transformed negative plugs", GH_ParamAccess.list);
            pManager.AddGeometryParameter("Sockets positive", "Sockets+", "Transformed positive sockets", GH_ParamAccess.list);
            pManager.AddGeometryParameter("Sockets negative", "Sockets-", "Transformed negative sockets", GH_ParamAccess.list);

            Params.Input[2].Optional = true;
            Params.Input[3].Optional = true;
            Params.Input[4].Optional = true;
            Params.Input[5].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGeometryParameter("Pieces with connectors", "PCs",
                "DataTree: each branch is one piece. Item 0 is the mesh, remaining items are its connectors.",
                GH_ParamAccess.tree);
            pManager.AddTextParameter("Piece names", "Ns", "Piece names in branch order", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            string jsonPath = "";
            if (!DA.GetData(0, ref jsonPath)) return;

            if (!File.Exists(jsonPath)) {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"JSON file not found: {jsonPath}");
                return;
            }

            string text;
            try {
                text = File.ReadAllText(jsonPath);
            } catch (Exception e) {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Could not read JSON: {e.Message}");
                return;
            }

            Data content;
            try {
                content = JsonSerializer.Deserialize<Data>(text);
            } catch (Exception e) {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Could not parse JSON: {e.Message}");
                return;
            }

            var meshGoos = new List<IGH_GeometricGoo>();
            if (!DA.GetDataList(1, meshGoos)) return;

            var plugsPos   = new List<IGH_GeometricGoo>(); DA.GetDataList(2, plugsPos);
            var plugsNeg   = new List<IGH_GeometricGoo>(); DA.GetDataList(3, plugsNeg);
            var socketsPos = new List<IGH_GeometricGoo>(); DA.GetDataList(4, socketsPos);
            var socketsNeg = new List<IGH_GeometricGoo>(); DA.GetDataList(5, socketsNeg);

            // build lookup: piece name -> connectors
            var connectorsByPiece = new Dictionary<string, List<IGH_GeometricGoo>>();
            foreach (var name in content.piece_names)
                connectorsByPiece[name] = new List<IGH_GeometricGoo>();

            void routeConnectors(List<IGH_GeometricGoo> goos) {
                foreach (var goo in goos) {
                    if (goo == null) continue;
                    var geo = goo.IsValid ? GH_Convert.ToGeometryBase(goo) : null;
                    if (geo == null) continue;
                    string pieceID = geo.GetUserString("pieceID");
                    if (pieceID != null && connectorsByPiece.ContainsKey(pieceID))
                        connectorsByPiece[pieceID].Add(goo);
                }
            }

            routeConnectors(plugsPos);
            routeConnectors(plugsNeg);
            routeConnectors(socketsPos);
            routeConnectors(socketsNeg);

            var tree = new GH_Structure<IGH_GeometricGoo>();
            var pieceNames = new List<string>();

            for (int i = 0; i < content.piece_names.Length; i++) {
                string pieceName = content.piece_names[i];
                var path = new GH_Path(i);

                // mesh at index i
                if (i < meshGoos.Count && meshGoos[i] != null)
                    tree.Append(meshGoos[i], path);
                else
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, $"No mesh at index {i} for piece {pieceName}");

                // connectors routed to this piece
                foreach (var goo in connectorsByPiece[pieceName])
                    tree.Append(goo, path);

                pieceNames.Add(pieceName);
            }

            DA.SetDataTree(0, tree);
            DA.SetDataList(1, pieceNames);
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("lego.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("9094345E-E59C-4E5D-B260-E24CFC11EFC3"); }
        }
    }
}