using System;
using System.Collections.Generic;
using System.Configuration;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text.Json;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;


///<summary>Component which creates and places connector geometry at singularities on a surface
///Placement is determined by br_piece_data.json which can be created using bertini_real write_piece_data()
///This component accepts 4 different connector geometry and requires at least 1 to run</summary>
namespace bertini_real
{
    public class SurfacePlaceComponents : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public SurfacePlaceComponents()
          : base("Surface Place Components", "SurfPlaceComps",
              "Read the specs for a surface from json, and place components at singularities, etc",
              "bertini_real", "Surface")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// They can be accessed in SolveInstance with DA.GetData()
        /// These will appear on the side of the component in which they appear
        /// can be accessed in the solve isntance in this order. 
        /// If you change their order here you MUST change the index used to access them in the SolveInstance
        /// DO NOT change the order of these once published/finalized
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            ///File to get data from. should be br_complete.json which is generated using write_piece() in python
            pManager.AddTextParameter("File Path", "F", "Path of json file with surface specs", GH_ParamAccess.tree);

            ///Size and location play to adjust connectors. No inputs required by user because it has a default value
            pManager.AddNumberParameter("Size", "S", "Scale factor for components", GH_ParamAccess.tree, 0.01);

            ///The connector Brep prefabs to place. At least one is required for the component to run, but it does not matter which one so all should be optional
            pManager.AddGeometryParameter("Plug Positive", "Plug+", "Plug positive geometry", GH_ParamAccess.tree);
            pManager.AddGeometryParameter("Plug Negative", "Plug-", "Plug negative geometry", GH_ParamAccess.tree);
            pManager.AddGeometryParameter("Socket Positive", "Socket+", "Socket positive geometry", GH_ParamAccess.tree);
            pManager.AddGeometryParameter("Socket Negative", "Socket-", "Socket negative geometry", GH_ParamAccess.tree);
            ///All the geometries should be optional. We check that there is at least 1 geo input in the SolveInstance
            Params.Input[1].Optional = true;
            Params.Input[2].Optional = true;
            Params.Input[3].Optional = true;
            Params.Input[4].Optional = true;
            Params.Input[5].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// can be set in the SolveInstance using DA.SetData()
        /// Appear on the side of the component in the order which they are listed
        /// Do NOT change their order once published/finalized
        /// if the order is changed you MUST UPDATE their index in the SolveInstance
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            ///Note: only the connector geometries with an inputted prefab should be sent to output
            ///Output a tree of the positive Brep connectors transformed to every singularity, one branch per piece
            pManager.AddGeometryParameter("Plugs positive", "Plugs+", "Pos geos transformed, per piece", GH_ParamAccess.tree);
            pManager.AddGeometryParameter("Plugs negative", "Plugs-", "Neg geos transformed, per piece", GH_ParamAccess.tree);
            ///Output a tree of the positive Brep connectors transformed to every singularity, one branch per piece
            pManager.AddGeometryParameter("Sockets positive", "Sockets+", "Pos geos transformed, per piece", GH_ParamAccess.tree);
            pManager.AddGeometryParameter("Sockets negative", "Sockets-", "Neg geos transformed, per piece", GH_ParamAccess.tree);
            ///Output a list of the piece filenames
            pManager.AddTextParameter("Piece filenames", "Fs", "Piece filenames", GH_ParamAccess.list);
            pManager.AddVectorParameter("Locations", "Locs", "Singularity locations as vectors", GH_ParamAccess.list);  // index 5
            pManager.AddVectorParameter("Directions", "Dirs", "Singularity directions as vectors", GH_ParamAccess.list); // index 6

            pManager.AddIntegerParameter("Sing Indices per Piece", "SingInds/Piece", "Singularity indices for each piece", GH_ParamAccess.tree); // index 7
            pManager.AddIntegerParameter("Sing Parities per Piece", "SingPars/Piece", "Singularity parities for each piece", GH_ParamAccess.tree); // index 8
            pManager.AddVectorParameter("Sing Directions per Piece", "SingDirs/Piece", "Singularity directions for each piece", GH_ParamAccess.tree); // index 9
            pManager.AddVectorParameter("Sing Locations per Piece", "SingLocs/Piece", "Singularity locations for each piece", GH_ParamAccess.tree); // index 10
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// GOAL: Take br_piece_data.json as produced from write_piece_data() in bertini_real
        /// Produce a pair of connectors at each singularity
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            /* Get all inputs
             * Create some empty variable
             * pass the data from the input parameters to the variables */

            ///empty variables
            Brep plugPos = null;
            Brep plugNeg = null;
            Brep socketPos = null;
            Brep socketNeg = null;
            Double size = 0.01;
            Point3d locationPlay = new Point3d(); // silviana sez: i don't think there's a setter for this.
            string jsonPath = "";

            List<String> piece_filenames = new List<String>();

            ///helper to extract a Brep from the first item of a geometry tree
            Brep BrepFromTree(GH_Structure<IGH_GeometricGoo> tree) {
                if (tree == null || tree.IsEmpty) return null;
                var goo = tree.get_FirstItem(true);
                if (goo == null) return null;
                // try direct cast first
                if (goo is GH_Brep ghBrep) return ghBrep.Value;
                // try casting via geometry base
                var geo = goo.IsValid ? GH_Convert.ToGeometryBase(goo) : null;
                if (geo is Brep brep) return brep;
                return null;
            }

            ///The parameters are stored in an array. 
            ///To set a variable to a parameter we need to reference the parameter by its index
            ///Do NOT want to change the order of these once published/finalized

            ///get file path from tree
            GH_Structure<GH_String> pathTree = new GH_Structure<GH_String>();
            DA.GetDataTree(0, out pathTree);
            if (pathTree.IsEmpty) return;
            jsonPath = pathTree.get_FirstItem(true).Value;

            ///get size from tree
            GH_Structure<GH_Number> sizeTree = new GH_Structure<GH_Number>();
            DA.GetDataTree(1, out sizeTree);
            if (!sizeTree.IsEmpty) size = sizeTree.get_FirstItem(true).Value;

            ///get geometry inputs from trees, take first item from each
            GH_Structure<IGH_GeometricGoo> plugPosTree   = new GH_Structure<IGH_GeometricGoo>();
            GH_Structure<IGH_GeometricGoo> plugNegTree   = new GH_Structure<IGH_GeometricGoo>();
            GH_Structure<IGH_GeometricGoo> socketPosTree = new GH_Structure<IGH_GeometricGoo>();
            GH_Structure<IGH_GeometricGoo> socketNegTree = new GH_Structure<IGH_GeometricGoo>();

            DA.GetDataTree(2, out plugPosTree);
            DA.GetDataTree(3, out plugNegTree);
            DA.GetDataTree(4, out socketPosTree);
            DA.GetDataTree(5, out socketNegTree);

            plugPos   = BrepFromTree(plugPosTree);
            plugNeg   = BrepFromTree(plugNegTree);
            socketPos = BrepFromTree(socketPosTree);
            socketNeg = BrepFromTree(socketNegTree);

            ///Error checking inputs. Including a RuntimeMessage in script will automatically generate an 'o' output on the component 

            ///The component should not run if there are no prefab geometries
            if ((plugNeg == null && plugPos != null) || (socketNeg == null && socketPos != null)) {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "Only positive geos inputted, ensure matching negative connectors are placed before combining with piece!"); 
            } //Remind the user if they only have positive geometries inputted that they will need negative geos if they want to combine with piece
            else if(plugNeg == null && plugPos == null && socketNeg == null && socketPos == null) {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "At least one of the plug/socket pos/neg geometries is invalid!");
                return;
            }

            /* Read and Parse the JSON File into a Data Object (defined in PlugParts.cs) */
            if (!File.Exists(jsonPath)) {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"JSON file not found: {jsonPath}");
                return;
            }

            string text;
            try {
                text = File.ReadAllText(jsonPath);
            } catch (Exception e) {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Could not read JSON file: {e.Message}");
                return;
            }

            ///this parses the JSON by key. The Data class must have properties the same name as the keys in the JSON file
            ///Should Eventually include some runtimeMessage error handling
            Data content;
            try
            {
                content = JsonSerializer.Deserialize<Data>(text);
            }
            catch (Exception e)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Could not parse JSON: {e.Message}");
                return;
            }

            ///JSON file is structured: 
            // {
            // "piece_names":['filename1.stl', 'filename2.stl', ...], 
            // "singularities_on_pieces"[[0,1],[2,3],[...], ...]
            // "sing_directions": [[vector for sing 0], [vec for sing 1], ... [vec for sing -1]],
            // "sing_locations": [[loc for sing 0], [loc for sing 1], ... [loc for sing -1]],
            // "parities": [[parities for singularity 0 on its pieces], [for sing 1], ..., [for sing -1]]
            // }
            // 
            ///where each property is a list of N lists, where N is the number of pieces. 
            ///Each list in a property corresponds to the property of the piece
            ///each piece is defined by the properties at the same index in each property list
            ///ex. piece 2 has piece_indices[1], singularities_on_piece[1], sing_directions[1], etc
            ///This is a silly way of using JSON because now we need to sort the JSON into each piece
            ///How we make the JSON in bertini_real write_piece_data really should be rewritten to organize by Sing or Piece where each Sing (or piece) has properties
            
            /* parse JSON data Piece objects */
            ///list for all the pieces
            
            // make the containers the the outputs
            ///one branch per piece, each branch contains the connectors for that piece
            GH_Structure<GH_Brep> plugs_pos_per_piece   = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Brep> plugs_neg_per_piece   = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Brep> sockets_pos_per_piece = new GH_Structure<GH_Brep>();
            GH_Structure<GH_Brep> sockets_neg_per_piece = new GH_Structure<GH_Brep>();

            GH_Structure<GH_Integer> sing_parities_per_piece = new GH_Structure<GH_Integer>();
            GH_Structure<GH_Integer> sing_indices_per_piece = new GH_Structure<GH_Integer>();
            GH_Structure<GH_Vector> sing_directions_per_piece = new GH_Structure<GH_Vector>();
            GH_Structure<GH_Vector> sing_locations_per_piece = new GH_Structure<GH_Vector>();

            GH_Structure<GH_Vector> all_sing_points = new GH_Structure<GH_Vector>(); // not per-piece.  just
            GH_Structure<GH_Vector> all_directions = new GH_Structure<GH_Vector>();

            List<PieceData> allPieces = new List<PieceData>();

            // we start with the per-piece structure
            // each piece is represented by a list of indices. the number of pieces = length of piece_indices
            for (int piece_index = 0; piece_index < content.piece_names.Length; piece_index++) {
                GH_Path path = new GH_Path(piece_index);

                plugs_pos_per_piece.EnsurePath(path);
                plugs_neg_per_piece.EnsurePath(path);
                sockets_pos_per_piece.EnsurePath(path);
                sockets_neg_per_piece.EnsurePath(path);

                string pieceName = content.piece_names[piece_index];

                piece_filenames.Add(pieceName);

                PieceData newPiece = new PieceData();
                newPiece.piece_name = pieceName; 
                // newPiece.indices = content.piece_indices[pieceName]; 
                newPiece.singsOnPiece = content.singularities_on_pieces[piece_index];
                
                //there are vectors for each sing on the piece, need to turn the vectors from vectors into lists
                //also append the vectors to the direction and location vector lists
                for (int j = 0; j < newPiece.singsOnPiece.Length; j++)
                {
                    int singIndex = newPiece.singsOnPiece[j];

                    sing_indices_per_piece.Append(new GH_Integer(singIndex), path);
                    sing_parities_per_piece.Append(new GH_Integer(content.parities[singIndex][piece_index]), path);

                    Vector3d direction = new Vector3d(
                        content.sing_directions[singIndex][0], 
                        content.sing_directions[singIndex][1], 
                        content.sing_directions[singIndex][2]);

                    Vector3d sing_point = new Vector3d(
                        content.sing_locations[singIndex][0], 
                        content.sing_locations[singIndex][1], 
                        content.sing_locations[singIndex][2]);
                    
                    sing_directions_per_piece.Append(new GH_Vector(direction), path);
                    sing_locations_per_piece.Append(new GH_Vector(sing_point), path);

                    /*Place correct connector on the piece at the singularity
                     * if the piece has positive polarity on the singularity, place a negative and positive Plug
                     * if it has negative polarity, add both socket pieces
                     * Only create new geometries if that geometry was given as an input */
                    


                    if (content.parities[singIndex][piece_index] == 1)
                    {
                        ///add the direction and location vectors for this plug to the list of all vectors and locations
                        ///these lists are now unused and can be deleted, but I am keeping them for debugging
                        /// 

                        // all_directions.Add(direction);
                        // all_sing_points.Add(sing_point);

                        if (plugNeg != null) { 
                            ///Create a new plug at this location and add it to the plug list, in the branch for this piece
                            plugs_neg_per_piece.Append(new GH_Brep(moveComponents(newPiece.piece_name, locationPlay, size, direction, sing_point, plugNeg)), path);
                        }
                        if (plugPos != null) {
                            ///Create a new plug at this location and add it to the plug list, in the branch for this piece
                            plugs_pos_per_piece.Append(new GH_Brep(moveComponents(newPiece.piece_name, locationPlay, size, direction, sing_point, plugPos)), path);
                        }
                    }

                    //add sockets if negative parity
                    else if (content.parities[singIndex][piece_index] == -1)
                    {
                        if (socketNeg != null)
                        {
                            ///Create a new socket at this location and add it to the socket list, in the branch for this piece
                            sockets_neg_per_piece.Append(new GH_Brep(moveComponents(newPiece.piece_name, locationPlay, size, direction, sing_point, socketNeg)), path);
                        }
                        if (socketPos != null)
                        {
                            ///Create a new socket at this location and add it to the socket list, in the branch for this piece
                            sockets_pos_per_piece.Append(new GH_Brep(moveComponents(newPiece.piece_name, locationPlay, size, direction, sing_point, socketPos)), path);
                        }                        
                    }
                }
            }

            for (int i = 0; i < content.sing_directions.Count(); i++)
            {
                GH_Path path = new GH_Path(0); // all in one branch since these are not per-piece

                Vector3d direction = new Vector3d(
                    content.sing_directions[i][0],
                    content.sing_directions[i][1],
                    content.sing_directions[i][2]);

                all_directions.Append(new GH_Vector(direction), path);

                Vector3d sing_point = new Vector3d(
                    content.sing_locations[i][0],
                    content.sing_locations[i][1],
                    content.sing_locations[i][2]);
                all_sing_points.Append(new GH_Vector(sing_point), path);
            }

            /* Set output data
             * 0 - Plugs positive tree (one branch per piece)
             * 1 - Plugs negative tree (one branch per piece)
             * 2 - Sockets positive tree (one branch per piece)
             * 3 - Sockets negative tree (one branch per piece)
             * 4 - Piece filenames list
             * 5 - Singularity locations list
             * 6 - Singularity directions list */
            DA.SetDataTree(0, plugs_pos_per_piece);
            DA.SetDataTree(1, plugs_neg_per_piece);
            DA.SetDataTree(2, sockets_pos_per_piece);
            DA.SetDataTree(3, sockets_neg_per_piece);
            DA.SetDataList(4, piece_filenames);

            // flattened things, not per-piece
            DA.SetDataList(5, all_sing_points);
            DA.SetDataList(6, all_directions);


            DA.SetDataTree(7, sing_indices_per_piece);
            DA.SetDataTree(8, sing_parities_per_piece);
            DA.SetDataTree(9, sing_directions_per_piece);
            DA.SetDataTree(10, sing_locations_per_piece);
        }

        /// <summary>
        /// Helper function which makes a new connector and places it at the singularity on the piece
        /// </summary>
        /// <param name="locationPlay">User input which changes the distance of the connector from the singularity</param>
        /// <param name="size">User input for scaling of the connector</param>
        /// <param name="direction">Direction vector points from the center of the piece to the singularity</param>
        /// <param name="location">Location of the singularity where the connector belongs</param>
        /// <param name="geo">Geometry of the connector prefab to be created</param>
        /// <returns>A new connector Brep at a singularity</returns>
        private Brep moveComponents(string piece_name, Point3d locationPlay, double size, Vector3d direction, Vector3d location, Brep geo) {
            ///create a new connector
            Brep newConnector = geo.DuplicateBrep(); 
            
            ///find our angles
            double phi = Math.Acos(direction[2] / direction.Length);
            double theta = Math.Atan2(direction[1], direction[0]);

            ///create some transformation matrices and then transform the connector
            var sf = Transform.Scale(Point3d.Origin + locationPlay, size);
            var rf = Transform.Rotation(phi, Vector3d.YAxis, Point3d.Origin); 
            
            newConnector.Transform(sf);
            newConnector.Transform(rf);
            
            rf = Transform.Unset; //clear the rotation matrix to be reused
            
            rf = Transform.Rotation(theta, Vector3d.ZAxis, Point3d.Origin);
            newConnector.Transform(rf);

            var xf = Transform.Translation(location);
            newConnector.Transform(xf);

            ///add user data so the piece the connector is attached to can later be identified
            newConnector.SetUserString("pieceID", piece_name);
            ///send back the connector
            return newConnector;
        }
        
        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("transform.png");

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("E636CDFC-C219-49A7-999A-06E91DE10B94"); }
        }
    }
}

///How to Importing STL
///String stlPiecePath = stlPath + "\\br_piece_smooth_"+piece.indices[0]+ "-" +piece.indices[1]+"-" +piece.indices[2]+"_solid.stl";
///DA.SetData(4, stlPiecePath);

///ActiveDoc.Import(stlPiecePath);