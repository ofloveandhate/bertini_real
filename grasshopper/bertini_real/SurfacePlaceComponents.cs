using System;
using System.Collections.Generic;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;


///<summary>Component which creates and places connector geometry at singularities on a surface.
///Placement comes from the singularity data emitted by "Surface Read GH JSON" (Sing Locations,
///Sing Directions, Sing Parities, Sing On Pieces) -- one JSON, one reader, no second file.
///Accepts up to 4 connector geometries and requires at least 1 to run.</summary>
namespace bertini_real
{
    public class SurfacePlaceComponents : GH_Component
    {
        public SurfacePlaceComponents()
          : base("Surface Place Components", "SurfPlaceComps",
              "Place plug/socket connectors at singularities, driven by Surface Read GH JSON's singularity outputs",
              "bertini_real", "Surface")
        {
        }

        /// <summary>
        /// Inputs.  The four singularity inputs wire directly from "Surface Read GH JSON".
        /// Do NOT reorder once published; the SolveInstance indexes by position.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Sing Locations", "SL", "Singularity locations (from Surface Read GH JSON)", GH_ParamAccess.list);
            pManager.AddVectorParameter("Sing Directions", "SD", "Singularity connector directions (from Surface Read GH JSON)", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Sing Parities", "SP", "Per singularity: parity (-1/0/1) on each piece (from Surface Read GH JSON)", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("Sing On Pieces", "SOP", "Per piece: indices of the singularities on it (from Surface Read GH JSON)", GH_ParamAccess.tree);

            // connector geometry is placed at true size -- orient + translate only, no scaling.
            pManager.AddGeometryParameter("Plug Positive", "Plug+", "Plug positive geometry", GH_ParamAccess.item);
            pManager.AddGeometryParameter("Plug Negative", "Plug-", "Plug negative geometry", GH_ParamAccess.item);
            pManager.AddGeometryParameter("Socket Positive", "Socket+", "Socket positive geometry", GH_ParamAccess.item);
            pManager.AddGeometryParameter("Socket Negative", "Socket-", "Socket negative geometry", GH_ParamAccess.item);

            Params.Input[4].Optional = true;  // geometries are individually optional; we require >=1
            Params.Input[5].Optional = true;
            Params.Input[6].Optional = true;
            Params.Input[7].Optional = true;
        }

        /// <summary>
        /// Outputs.  One branch per piece for the connector trees and per-piece diagnostics.
        /// Do NOT reorder once published.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGeometryParameter("Plugs positive", "Plugs+", "Positive plugs transformed, per piece", GH_ParamAccess.tree);
            pManager.AddGeometryParameter("Plugs negative", "Plugs-", "Negative plugs transformed, per piece", GH_ParamAccess.tree);
            pManager.AddGeometryParameter("Sockets positive", "Sockets+", "Positive sockets transformed, per piece", GH_ParamAccess.tree);
            pManager.AddGeometryParameter("Sockets negative", "Sockets-", "Negative sockets transformed, per piece", GH_ParamAccess.tree);

            pManager.AddIntegerParameter("Sing Indices per Piece", "SingInds/Piece", "Singularity indices for each piece", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("Sing Parities per Piece", "SingPars/Piece", "Singularity parities for each piece", GH_ParamAccess.tree);
            pManager.AddVectorParameter("Sing Directions per Piece", "SingDirs/Piece", "Singularity directions for each piece", GH_ParamAccess.tree);
            pManager.AddPointParameter("Sing Locations per Piece", "SingLocs/Piece", "Singularity locations for each piece", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            var locations = new List<Point3d>();
            var directions = new List<Vector3d>();
            GH_Structure<GH_Integer> parities;
            GH_Structure<GH_Integer> onPieces;

            DA.GetDataList(0, locations);
            DA.GetDataList(1, directions);
            if (!DA.GetDataTree(2, out parities)) return;
            if (!DA.GetDataTree(3, out onPieces)) return;

            Brep plugPos = BrepFromInput(DA, 4);
            Brep plugNeg = BrepFromInput(DA, 5);
            Brep socketPos = BrepFromInput(DA, 6);
            Brep socketNeg = BrepFromInput(DA, 7);

            if (plugPos == null && plugNeg == null && socketPos == null && socketNeg == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "At least one plug/socket geometry is required.");
                return;
            }
            if ((plugPos != null && plugNeg == null) || (socketPos != null && socketNeg == null))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Remark,
                    "Only positive geometry supplied for a connector; add the matching negative before combining with the piece.");
            }

            var plugsPos = new GH_Structure<GH_Brep>();
            var plugsNeg = new GH_Structure<GH_Brep>();
            var socketsPos = new GH_Structure<GH_Brep>();
            var socketsNeg = new GH_Structure<GH_Brep>();
            var indicesPerPiece = new GH_Structure<GH_Integer>();
            var paritiesPerPiece = new GH_Structure<GH_Integer>();
            var dirsPerPiece = new GH_Structure<GH_Vector>();
            var locsPerPiece = new GH_Structure<GH_Vector>();

            // one branch per piece, from the Sing On Pieces tree
            for (int b = 0; b < onPieces.PathCount; b++)
            {
                GH_Path piecePath = onPieces.get_Path(b);
                int pieceIndex = piecePath.Indices.Length > 0 ? piecePath.Indices[piecePath.Indices.Length - 1] : b;

                plugsPos.EnsurePath(piecePath);
                plugsNeg.EnsurePath(piecePath);
                socketsPos.EnsurePath(piecePath);
                socketsNeg.EnsurePath(piecePath);

                foreach (var goo in onPieces.get_Branch(piecePath))
                {
                    if (!(goo is GH_Integer gi)) continue;
                    int singIndex = gi.Value;

                    if (singIndex < 0 || singIndex >= locations.Count || singIndex >= directions.Count)
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, $"Singularity index {singIndex} out of range.");
                        continue;
                    }

                    int parity = ParityFor(parities, singIndex, pieceIndex);
                    Vector3d direction = directions[singIndex];
                    Point3d location = locations[singIndex];

                    indicesPerPiece.Append(new GH_Integer(singIndex), piecePath);
                    paritiesPerPiece.Append(new GH_Integer(parity), piecePath);
                    dirsPerPiece.Append(new GH_Vector(direction), piecePath);
                    locsPerPiece.Append(new GH_Vector(new Vector3d(location)), piecePath);

                    // parity +1 -> this piece gets the plug; -1 -> it gets the socket
                    if (parity == 1)
                    {
                        if (plugNeg != null)
                            plugsNeg.Append(new GH_Brep(moveComponents(pieceIndex.ToString(), direction, new Vector3d(location), plugNeg)), piecePath);
                        if (plugPos != null)
                            plugsPos.Append(new GH_Brep(moveComponents(pieceIndex.ToString(), direction, new Vector3d(location), plugPos)), piecePath);
                    }
                    else if (parity == -1)
                    {
                        if (socketNeg != null)
                            socketsNeg.Append(new GH_Brep(moveComponents(pieceIndex.ToString(), direction, new Vector3d(location), socketNeg)), piecePath);
                        if (socketPos != null)
                            socketsPos.Append(new GH_Brep(moveComponents(pieceIndex.ToString(), direction, new Vector3d(location), socketPos)), piecePath);
                    }
                }
            }

            DA.SetDataTree(0, plugsPos);
            DA.SetDataTree(1, plugsNeg);
            DA.SetDataTree(2, socketsPos);
            DA.SetDataTree(3, socketsNeg);
            DA.SetDataTree(4, indicesPerPiece);
            DA.SetDataTree(5, paritiesPerPiece);
            DA.SetDataTree(6, dirsPerPiece);
            DA.SetDataTree(7, locsPerPiece);
        }

        /// <summary>Parity (-1/0/1) of a singularity on a given piece, from the per-singularity parity tree.</summary>
        private static int ParityFor(GH_Structure<GH_Integer> parities, int singIndex, int pieceIndex)
        {
            GH_Path path = new GH_Path(singIndex);
            if (!parities.PathExists(path)) return 0;
            var branch = parities.get_Branch(path);
            if (pieceIndex < 0 || pieceIndex >= branch.Count) return 0;
            return (branch[pieceIndex] as GH_Integer)?.Value ?? 0;
        }

        /// <summary>Pull a Brep from an item geometry input, or null if absent/unconvertible.</summary>
        private static Brep BrepFromInput(IGH_DataAccess DA, int index)
        {
            IGH_GeometricGoo goo = null;
            if (!DA.GetData(index, ref goo) || goo == null) return null;
            var geo = GH_Convert.ToGeometryBase(goo);
            return geo as Brep;
        }

        /// <summary>
        /// Make a new connector and place it at a singularity: orient to the direction
        /// (phi about Y, theta about Z) and translate to the location.  The connector keeps its
        /// own (true) size -- no scaling here; size it upstream.
        /// </summary>
        private Brep moveComponents(string piece_name, Vector3d direction, Vector3d location, Brep geo)
        {
            Brep newConnector = geo.DuplicateBrep();

            double phi = Math.Acos(direction[2] / direction.Length);
            double theta = Math.Atan2(direction[1], direction[0]);

            var rf = Transform.Rotation(phi, Vector3d.YAxis, Point3d.Origin);
            newConnector.Transform(rf);

            rf = Transform.Rotation(theta, Vector3d.ZAxis, Point3d.Origin);
            newConnector.Transform(rf);

            var xf = Transform.Translation(location);
            newConnector.Transform(xf);

            newConnector.SetUserString("pieceID", piece_name);
            return newConnector;
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("transform.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("E636CDFC-C219-49A7-999A-06E91DE10B94"); }
        }
    }
}
