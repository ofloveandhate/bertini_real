using System;
using System.Collections;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Collections;
using Rhino.Geometry;
using Rhino.PlugIns;
using static Rhino.Render.TextureGraphInfo;

namespace MyPlugin.UtilityComps
{
    public class Connectors : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public Connectors()
          : base("Connnectors ", "Conn",
              "All of the connetor prefabs for for nodes",
              "MyPlugin", "Utilites")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            ///PLUG PARAMS
            pManager.AddNumberParameter("Wire Hole Diameter", "D", "diameter of hole required for a wire", GH_ParamAccess.item, 7.0);
            pManager.AddNumberParameter("Plug Length Factor", "PF", "multiplier for length of plug neg", GH_ParamAccess.item, 3.40);
            pManager.AddNumberParameter("Wedge Height", "WH", "height of tab wedge", GH_ParamAccess.item, 4.0);
            pManager.AddNumberParameter("Wedge Length", "WL", "length of tab wedge", GH_ParamAccess.item, 2.0);
            pManager.AddNumberParameter("Wedge Width", "WW", "width of tab wedge", GH_ParamAccess.item, 2.0);
            pManager.AddNumberParameter("Cutout Thickness", "CT", "space between tab and center of plug", GH_ParamAccess.item, 1.40);
            pManager.AddNumberParameter("Cutout Depth", "CD", "how far down do the tabs go", GH_ParamAccess.item, 10.0);
            pManager.AddNumberParameter("Tab Thickness", "TT", "thickness of both tabs", GH_ParamAccess.item, 2.76);
            pManager.AddNumberParameter("Taper Factor", "TF", "how much smaller the tapered radius is to the plug radius", GH_ParamAccess.item, 0.70);
            pManager.AddNumberParameter("Taper Length", "TL", "distance between the smaller radius and main body cyllinder", GH_ParamAccess.item, 5.0);
            ///SOCKET PARAMS
            pManager.AddNumberParameter("Socket Diameter", "SD", "Diameter of socket", GH_ParamAccess.item, 15.00);
            pManager.AddNumberParameter("Socket Wall Thickness", "ST", "thickness from the outer diameter of the socket inward", GH_ParamAccess.item, 2.00);
            pManager.AddNumberParameter("Socket Length", "SL", "height of the socket", GH_ParamAccess.item, 7.0);
            pManager.AddNumberParameter("Socket Length Factor", "SF", "height of the socket", GH_ParamAccess.item, 3.0);
            ///SHARED PARAMS
            pManager.AddNumberParameter("Connection Play", "CP", "an adjustable amount to tweak the snuggness of the plugs fit into the socket", GH_ParamAccess.item, 0.20);
            pManager.AddNumberParameter("Body Overlap", "BO", "how much the plug and socket overlap", GH_ParamAccess.item, 7.0);
            pManager.AddNumberParameter("Length Overage", "LO", "how much excess hangover between the plug and socket", GH_ParamAccess.item, 7.0);
            pManager.AddNumberParameter("eps", "E", "an adjusment value default 0.01", GH_ParamAccess.item, 0.1);
            pManager.AddPlaneParameter("Base", "B", "base xy plane to build upon", GH_ParamAccess.item, Plane.WorldXY);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGeometryParameter("Negative Plug", "NP", "cylinder brep to be subtracted from the surface", GH_ParamAccess.item);
            pManager.AddGeometryParameter("Positive Plug", "PP", "cylinder brep to be subtracted from the surface", GH_ParamAccess.item);
            pManager.AddGeometryParameter("Negative Socket", "NS", "cylinder brep to be subtracted from the surface", GH_ParamAccess.item);
            pManager.AddGeometryParameter("Positive Socket", "PS", "cylinder brep to be subtracted from the surface", GH_ParamAccess.item);
            pManager.AddGeometryParameter("Debug", "NP", "cylinder brep to be subtracted from the surface", GH_ParamAccess.item);
        }

        
        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double wireHoleDia = 0;
            double plugFactor = 0;
            double wedgeHeight = 0;
            double wedgeLength = 0;
            double wedgeWidth = 0;
            double cutoutThickness = 0;
            double cutoutDepth = 0;
            double tabThickness = 0;
            double taperFactor = 0;
            double taperLength = 0;

            double socketDiameter = 0;
            double socketWallThickeness = 0;
            double socketLength = 0;
            double socketFactor = 0;

            double connectionPlay = 0;
            double bodyOverlap = 0;
            double lengthOverage = 0;
            double eps = 0;
            Plane b = Plane.Unset;

            if (!DA.GetData(0, ref wireHoleDia)) return;
            if (!DA.GetData(1, ref plugFactor)) return;
            if (!DA.GetData(2, ref wedgeHeight)) return;
            if (!DA.GetData(3, ref wedgeLength)) return;
            if (!DA.GetData(4, ref wedgeWidth)) return;
            if (!DA.GetData(5, ref cutoutThickness)) return;
            if (!DA.GetData(6, ref cutoutDepth)) return;
            if (!DA.GetData(7, ref tabThickness)) return;
            if (!DA.GetData(8, ref taperFactor)) return;
            if (!DA.GetData(9, ref taperLength)) return;
            if (!DA.GetData(10, ref socketDiameter)) return;
            if (!DA.GetData(11, ref socketWallThickeness)) return;
            if (!DA.GetData(12, ref socketLength)) return;
            if (!DA.GetData(13, ref socketFactor)) return;
            if (!DA.GetData(14, ref connectionPlay)) return;
            if (!DA.GetData(15, ref bodyOverlap)) return;
            if (!DA.GetData(16, ref lengthOverage)) return;
            if (!DA.GetData(17, ref eps)) return;
            if (!DA.GetData(18, ref b)) return;

            /**Build Positive plug*/

            double plugR = (socketDiameter / 2) - socketWallThickeness- connectionPlay; //radius of the plug shuld fit snug into the socket
            double taperedR = taperFactor * plugR;
            //main body height fits entirely into the socket plus some overlap and excess length, minus the taper length
            double bodyHeight = socketLength + bodyOverlap - taperLength + lengthOverage;

            /**Send data to objects for creation*/
            PlugBody plugBody = new PlugBody(taperedR, taperLength, plugR, bodyHeight); //used to create the plug body
            PlugTabs plugTabs = new PlugTabs(wedgeHeight, wedgeLength, wedgeWidth, cutoutThickness, cutoutDepth, tabThickness); //used to create cuoutBoxes and Wedges

            /** create the overall plug body*/
            Brep plugBodyBrep = plugBody.plugBodyGeo();
            //DA.SetData(0, plugBodyBrep); //send the plug body geometry to "plugBody" output for debugging

            /**Create and place cutout boxes to create the tabs
             * The difference of the boxes will be taken from the plugbody
             * Cutout boxes should be equi distance from the center of the plug
             */
            // var mf = Transform.Mirror(Plane.WorldZX); //create mirror transformation vector matrix
            Brep[] cutouts = new Brep[2] { plugTabs.cutoutBox(plugR, b, eps), plugTabs.cutoutBox(plugR, b, eps) }; //create two new cutout boxes
            plugBodyBrep = cutoutBoxes(cutouts, plugBodyBrep);


            /**Create and place wedges onto outside of plug body
             * Vertical Face should be parrallel to cutout box long face*/
            Brep[] wedges = new Brep[2] { plugTabs.rawWedge(), plugTabs.rawWedge() }; //Create 2 wedges
            ;
            plugBodyBrep = addWedges(wedges, plugR, tabThickness, taperLength, plugBodyBrep);
            DA.SetData(1, plugBodyBrep);


            /**Negative Socket
            */
            //Point3d origin = rg.Point3d(0, 0, 0)
            //yaxis = rg.Vector3d(0, 1, 0)
            double scaledSocketLength = socketFactor * socketLength + eps;
            // length of cyllinder is the length of the socket plus some excess, scaled up the socket_neg_factor
            //create the the socket as a cyllinder from a base circle that has a radius of the inner socket 
            Brep negativeSocket = new Cylinder(new Circle(b, (socketDiameter / 2) - socketWallThickeness), scaledSocketLength).ToBrep(true, true);
            //rotate and move the cyllinder into correct position
            var xf = Transform.Translation(0, 0, -scaledSocketLength / 2);
            negativeSocket.Transform(xf);
            var rf = Transform.Rotation(Math.PI, Vector3d.YAxis, Point3d.Origin);
            negativeSocket.Transform(rf);
            DA.SetData(2, negativeSocket);

            //unset transformation vectors
            xf=Transform.Unset;
            rf = Transform.Unset;

            /**Positive Socket*/
            Brep innerCylinder = new Cylinder(new Circle(b, (socketDiameter / 2) - socketWallThickeness), socketLength + 2 * eps).ToBrep(true, true);
            xf = Transform.Translation(0, 0, -eps);
            Brep outerCylinder = new Cylinder(new Circle(b, socketDiameter / 2), socketLength).ToBrep(true, true);

            Brep[] socket = Brep.CreateBooleanDifference(outerCylinder, innerCylinder, 0.01);
            rf = Transform.Rotation(Math.PI, Vector3d.YAxis, Point3d.Origin);
            socket[0].Transform(rf);
            DA.SetData(3, socket[0]);

            //unset transformation vectors
            xf = Transform.Unset;
            rf = Transform.Unset;

            /** Negative Plug*/
            double cuttingR = wireHoleDia / 2;
            double scaledPlugLength = plugFactor * (socketLength + lengthOverage + bodyOverlap);
            //create a cyllinder and move it down the z axis to center it on the origin
            Brep negativePlug = new Cylinder(new Circle(b, cuttingR), scaledPlugLength).ToBrep(true, true);
            xf = Transform.Translation(0, 0, -eps - (scaledPlugLength / 2));
            negativePlug.Transform(xf);
            DA.SetData(0, negativePlug);
            DA.SetData(4,negativePlug);
        }

        private Brep addWedges(Brep[] wedges, double plugR, double tabThickness, double plugTaperLength, Brep plugBodyBrep)
        {
            /**Create and place wedges onto outside of plug body
             * Vertical Face should be parrallel to cutout box long face*/
            var xf = Transform.Translation(0, plugR - tabThickness, plugTaperLength); //create transformation matrix for both wedges
            wedges[0].Transform(xf); //move wedges into place

            wedges[1].Transform(xf);
            wedges[1].Transform(Transform.Mirror(Plane.WorldZX)); //Flip one wedge to the other side of the plug

            //Brep[] thingsToUnion = new Brep[2] { plugBodyBrep, wedges[0] };
            //plugBodyBrep = Brep.CreateBooleanUnion(thingsToUnion, 0.01)[0];

            wedges[1].Faces.SplitKinkyFaces(RhinoMath.DefaultAngleTolerance, true);
            if (BrepSolidOrientation.Inward == wedges[1].SolidOrientation)
            {
                wedges[1].Flip();
            }
            Brep[] thingsToUnion = new Brep[3] { plugBodyBrep, wedges[0], wedges[1] };

            //thingsToUnion[1] = wedges[1];
            plugBodyBrep = Brep.CreateBooleanUnion(thingsToUnion, 0.01)[0];
            return plugBodyBrep;
        }
        private Brep cutoutBoxes(Brep[] cutouts, Brep plugBodyBrep)
        {
            /**Create and place cutout boxes to create the tabs
                 * The difference of the boxes will be taken from the plugbody
                 * Cutout boxes should be equi distance from the center of the plug
                 */
            var mf = Transform.Mirror(Plane.WorldZX); //create mirror transformation vector matrix
                                                      //Brep[] cutouts = new Brep[2] { plugTabs.cutoutBox(plugR, b, eps), plugTabs.cutoutBox(plugR, b, eps) }; //create two new cutout boxes
                                                      //DA.SetData(3, cutouts[0]); //send the cutout geometry to outparamter "cutout1" for debugging
            cutouts[1].Transform(mf); //mirror the second cutout to the other side of the plug body
                                      //DA.SetData(4, cutouts[1]); //send the cutout geometry to outparamter "cutout2" for debugging

            /**Cutout the cutout boxes from the plug body to create tabs*/
            plugBodyBrep = Brep.CreateBooleanDifference(plugBodyBrep, cutouts[0], 0.01)[0]; //plugBody-cutout[0]
                                                                                            //DA.SetData(1, plugBodyBrep); //Send plugBody with one cutout to "PBOneCut" for debugging 
            ///Because of the mirror, the polyline orientation of cutout[1] can become flipped from counter-clockwise to clockwise
            ///Boolean operations can become flipped when their orientation is flipped
            ///Read more: https://discourse.mcneel.com/t/rhino-geometry-brep-createbooleandifference-not-creating-difference/116008
            ///Solution: https://discourse.mcneel.com/t/rhino-geometry-brep-createbooleandifference-not-creating-difference/116008
            ///Without these 3 lines, the difference results in cutout-plugBody, not plugbody-cutout
            cutouts[1].Faces.SplitKinkyFaces(RhinoMath.DefaultAngleTolerance, true);
            if (BrepSolidOrientation.Inward == cutouts[1].SolidOrientation)
            {
                cutouts[1].Flip();
            }
            plugBodyBrep = Brep.CreateBooleanDifference(plugBodyBrep, cutouts[1], 0.01)[0]; //Cutout the second cutout box. Need Split KinkyFaces and Flip Orientation Lines! 
            return plugBodyBrep; //DA.SetData(2, plugBodyBrep); //send the plug body with both cutouts to "PBTwoCut" for debugging

        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return Properties.Resources.lego;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("9b510117-7c7b-4848-8dcc-e694003de5a7"); }
        }
    }
}