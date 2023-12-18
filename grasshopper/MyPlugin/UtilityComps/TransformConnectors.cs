using System;
using System.Collections.Generic;
using System.Configuration;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text.Json;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

///<summary>Component which creates and places connector geometry at singularities on a surface
///Placement is determined by br_piece_data.json which can be created using bertini_real write_piece_data()
///This component accepts 4 different connector geometry and requires at least 1 to run</summary>
namespace MyPlugin.UtilityComps
{
    public class TransformConnectors : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public TransformConnectors()
          : base("Transform Connectors", "Transform",
              "Description",
              "MyPlugin", "Utilites")
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
            pManager.AddTextParameter("File Path", "F", "Path of br complete file", GH_ParamAccess.item);
            ///Size and location play to adjust connectors. No inputs required by user because it has a default value
            pManager.AddNumberParameter("Size", "S", "Scale factor for components", GH_ParamAccess.item, 0.01);
            ///If no Location play is given default to no location Play (origin)
            pManager.AddPointParameter("Location Play", "L", "direction to adjust the connectors uniformly", GH_ParamAccess.item, Point3d.Origin);
            ///The connector Brep prefabs place. At least one is required for the component to run, but it does not matter which one so all should be optional
            pManager.AddGeometryParameter("Neg Plug", "NP", "Connector Prefab to place and tansform", GH_ParamAccess.item);
            pManager.AddGeometryParameter("Pos Plug", "PP", "Connector Prefab to place and tansform", GH_ParamAccess.item);
            pManager.AddGeometryParameter("Neg Socket", "NS", "Connector Prefab to place and tansform", GH_ParamAccess.item);
            pManager.AddGeometryParameter("Pos Socket", "PS", "Connector Prefab to place and tansform", GH_ParamAccess.item);
            ///All the geometries should be optional. We check that there is at least 1 geo input in the SolveInstance
            Params.Input[3].Optional = true;
            Params.Input[4].Optional = true;
            Params.Input[5].Optional = true;
            Params.Input[6].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// can be set in the SolveInstance using DA.SetData()
        /// Appear on the side of the component in the order which they are listed
        /// Do NOT change their order once published/finalized
        /// if the order is chang you MUST UPDATE their index in the SolveInstance
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            ///Note: only the connector geometries with an inputted prefab should be sent to output
            ///Output a list of the negative Brep connectors transformed to every singularity
            pManager.AddGeometryParameter("Neg Plugs", "NP", "negative transformed geos", GH_ParamAccess.list);
            ///Output a list of the positive Brep connectors transformed to every singularity
            pManager.AddGeometryParameter("Pos Plugs", "PP", "Pos geos transofrmed", GH_ParamAccess.list);
            pManager.AddGeometryParameter("Neg Socket", "NS", "negative transformed geos", GH_ParamAccess.list);
            ///Output a list of the positive Brep connectors transformed to every singularity
            pManager.AddGeometryParameter("Pos Socket", "PS", "Pos geos transofrmed", GH_ParamAccess.list);
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
            Brep negPlug = new Brep();
            Brep posPlug = new Brep();
            Brep negSocket = new Brep();
            Brep posSocket= new Brep();
            Double size = 0;
            Point3d locationPlay = new Point3d();
            List<Vector3d> locVectors = new List<Vector3d>();
            List<Vector3d> dirVectors = new List<Vector3d>();
            string jsonPath = "";
            
            List<Brep> transformedNegPlugs = new List<Brep>();
            List<Brep> transformedPosPlugs = new List<Brep>();
            List<Brep> transformedNegSockets = new List<Brep>();
            List<Brep> transformedPosSockets = new List<Brep>();

            ///The parameters are stored an array. 
            ///To set a variable to a parameter we need to reference the parameter by its index
            ///Do NOT want to change the order of these once published/finalized 
            if (!DA.GetData(0, ref jsonPath)) return;
            if (!DA.GetData(1, ref size)) return;
            if (!DA.GetData(2, ref locationPlay)) return;
            DA.GetData(3, ref negPlug);
            DA.GetData(4, ref posPlug);
            DA.GetData(5, ref negSocket);
            DA.GetData(6, ref posSocket);

            ///Error checking inputs. Including a RuntimeMessage in script will automaticall generate an 'o' output on the component 
            ///The component should not run in there are no prefab geometries
            if ((!negPlug.IsValid && posPlug.IsValid) || (!negSocket.IsValid && posSocket.IsValid)) {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "Only positive geos inputted, ensure matching negative connectors are placed before combining with piece!"); 
            } //Remind the user if they only have positive geometries inputted that they will need negative geos if they want to combine with piece
            else if(!negPlug.IsValid && !posPlug.IsValid && !negSocket.IsValid && !posSocket.IsValid) {
                this.AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "No input geos! Include at least 1 geo to transform");
                return;
            }

            /* Read and Parse the JSON File into a Data Object (defined in PlugParts.cs) */
            string text = File.ReadAllText(jsonPath);
            ///this parses the JSON by key. The Data class must have properties the same name as the keys in the JSON file
            ///Should Eventually include some runtimeMessage error handeling
            var content = JsonSerializer.Deserialize<Data>(text); 
            ///JSON file is structured: {piece_indices:[[i1,i2,i3,i4],[i1,i2,i3,i4],[...]], singularities_on_pieces[[s1,s2,s3],[s1,s2,s3],[...], ...}
            ///where each property is a list of N lists, where N is the number of pieces. 
            ///Each list  in a property corresponds to the property of the piece
            ///each piece is defined by the properties at the same index in each property list
            ///ex. piece 2 has piece_indices[1], singularities_on_piece[1], sing_directions[1], etc
            ///This is a silly way of using JSON because now we need to sort the JSON into each piece
            ///How we make the JSON in bertini_real write_piece_data really should be rewritten to organize by Sing or Piece where each Sing (or piece) has properties
            
            /* parse JSON data Piece objects */
            ///list for all the pieces
            List<PieceData> allPieces = new List<PieceData>();
            ///each piece is represented by a list of indices. the number of peices = length of piece_indices
            for (int pieceIndex= 0; pieceIndex < content.piece_indices.Length; pieceIndex++) {
                //⚠️I would like to try just pass content to PieceData and have it do the work for me!
                PieceData newPiece = new PieceData();
                newPiece.piece_index = pieceIndex; 
                newPiece.indices = content.piece_indices[pieceIndex]; 
                newPiece.singsOnPiece = content.singularities_on_pieces[pieceIndex];
                
                //there are vectors for each sing on the piece, need to turn the vectors from vectors into lists
                //also append the vectors to the direction and location vector lists
                for (int j = 0; j <newPiece.singsOnPiece.Length ; j++)
                {
                    int singIndex = newPiece.singsOnPiece[j];

                    ///these vector lists are now defunct but kept for setimental and debugging
                    Vector3d dirVect = new Vector3d(content.sing_directions[singIndex][0], content.sing_directions[singIndex][1], content.sing_directions[singIndex][2]);
                    Vector3d locVect = new Vector3d(content.sing_locations[singIndex][0], content.sing_locations[singIndex][1], content.sing_locations[singIndex][2]);
                    
                    /*Place correct connector on the piece at the singularity
                     * if the piece has positive polarity on the singularity, place a negative and positive Plug
                     * if it has negative polarity, add both socket pieces
                     * Only create new geometries if that geometry was given as an input */
                    if (content.parities[singIndex][pieceIndex] == 1)
                    {
                        ///add the direction and location vectors for this plug to the list of all vectors and locations
                        ///these lists are now unused and can be deleted, but I am keeping them for debugging
                        dirVectors.Add(dirVect);
                        locVectors.Add(locVect);

                        if (negPlug.IsValid) { 
                            ///Create a new plug at this location and add it to the plug list
                            transformedNegPlugs.Add(moveComponents(newPiece.indices,locationPlay,size,dirVect,locVect,negPlug));
                        }
                        if (posPlug.IsValid) {
                            ///Create a new plug at this location and add it to the plug list
                            transformedPosPlugs.Add(moveComponents(newPiece.indices, locationPlay, size, dirVect, locVect, posPlug));
                        }
                    }

                    //add sockets if negative parity
                    else if (content.parities[singIndex][pieceIndex] == -1)
                    {
                        ///add the direction and location vectors for this socket to the list of all vectors and locations
                        ///these lists are now unused and can be deleted, but I am keeping them for debugging
                        dirVectors.Add(dirVect);
                        locVectors.Add(locVect);
                        if (negSocket.IsValid)
                        {
                            ///Create a new plug at this location and add it to the plug list
                            transformedNegSockets.Add(moveComponents(newPiece.indices, locationPlay, size, dirVect, locVect, negSocket));
                        }
                        if (posSocket.IsValid)
                        {
                            ///Create a new plug at this location and add it to the plug list
                            transformedPosSockets.Add(moveComponents(newPiece.indices, locationPlay, size, dirVect, locVect, posSocket));
                        }                        
                    }
                }
            }

            
            
            /* Set output data
             * Set to the list of Geos
             * 0 - Out must be text
             * 1 - negComponents List
             * 2 - posCompoents List */
            DA.SetDataList(0, transformedNegPlugs);
            DA.SetDataList(1, transformedPosPlugs);
            DA.SetDataList(2, transformedNegSockets);
            DA.SetDataList(3, transformedPosSockets);
        }

        /// <summary>
        /// Helper function which makes a new connector and places it at the singularity on the piece*/
        /// </summary>
        /// <param name="locationPlay">User input which changes the distance of the connector from the singularity</param>
        /// <param name="size">User input for scaling of the connector</param>
        /// <param name="direction">Direction vector points from the center of the piece to the singularity</param>
        /// <param name="location">Location of the singularity where the connector belongs</param>
        /// <param name="geo">Geometry of the connector prefab to be created</param>
        /// <returns>A new connector Brep at a singularity</returns>
        private Brep moveComponents(int[] pieceIndices, Point3d locationPlay, double size, Vector3d direction, Vector3d location, Brep geo) {
            ///create a new connector
            Brep newConnector = geo.DuplicateBrep(); 
            
            ///find our angles
            double phi = Math.Acos(direction[2] / direction.Length);
            double theta = Math.Atan2(direction[1], direction[0]);

            ///create some transformation matricies and then tranform the connector
            var sf = Transform.Scale(Point3d.Origin + locationPlay, size);
            var rf = Transform.Rotation(phi, Vector3d.YAxis, Point3d.Origin); 
            
            newConnector.Transform(sf);
            newConnector.Transform(rf);
            
            rf = Transform.Unset; //clear the rotation matrix to be reused
            
            rf = Transform.Rotation(theta, Vector3d.ZAxis, Point3d.Origin);
            newConnector.Transform(rf);

            var xf = Transform.Translation(location);
            newConnector.Transform(xf);

            ///add user data so the piece the the connector is attached to can later be identified
            newConnector.SetUserString("pieceID", pieceIndices.ToString());
            ///send back the connector
            return newConnector;

        }
        
        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                return Properties.Resources.transform_icon;
            }
        }

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