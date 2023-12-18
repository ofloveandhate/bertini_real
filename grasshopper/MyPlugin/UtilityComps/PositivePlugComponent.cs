using System;
using System.Collections;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino;
using Rhino.Collections;
using Rhino.Geometry;
using Rhino.Geometry.Collections;
using Rhino.PlugIns;
using Rhino.Render.ChangeQueue;
using Rhino.Render.PostEffects;

//public class PlugBody
//{
//    public double taperedR { get; set; }
//    public double plugTaperLength { get; set; }
//    public double plugR { get; set; }
//    public double bodyHeight { get; set; }
//    public PlugBody(double taperedR, double plugTaperLength, double plugR, double bodyHeight)
//    {
//        this.taperedR = taperedR;
//        this.plugTaperLength = plugTaperLength;
//        this.plugR = plugR;
//        this.bodyHeight = bodyHeight;
//    }
//    public Brep plugBodyGeo()
//    {
//        Point3d p0 = new Point3d(0, 0, 0);
//        Point3d p1 = new Point3d(taperedR, 0, 0); //origin to taper radius
//        Point3d p2 = new Point3d(plugR, 0, plugTaperLength); //transition from taper to body
//        Point3d p3 = new Point3d(plugR, 0, plugTaperLength + bodyHeight); //to top of plug
//        Point3d p4 = new Point3d(0, 0, plugTaperLength + bodyHeight); //top on axis

//        Polyline polyline = new Polyline(new Point3dList(p0, p1, p2, p3));
//        Line axis = new Line(p0, p4);
//        Brep geo = RevSurface.Create(polyline, axis).ToBrep();
//        return geo.CapPlanarHoles(0.001);
//    }
//}

//public class PlugTabs
//{
//    public double wedgeHeight { get; set; }
//    public double wedgeLength { get; set; }
//    public double wedgeWidth { get; set; }
//    public double tabCutoutThickness { get; set; }
//    public double tabCutoutDepth { get; set; }
//    public double tabThickness { get; set; }

//    public PlugTabs(double wedgeHeight, double wedgeLength, double wedgeWidth, double tabCutoutThickness, double tabCutoutDepth, double tabThickness)
//    {
//        this.wedgeHeight = wedgeHeight;
//        this.wedgeLength = wedgeLength;
//        this.wedgeWidth = wedgeWidth;
//        this.tabCutoutThickness = tabCutoutThickness;
//        this.tabCutoutDepth = tabCutoutDepth;
//        this.tabThickness = tabThickness;
//    }
//    public Brep cutoutBox(double plugR, Plane b, double eps)
//    {
//        double ell = plugR - tabThickness - tabCutoutThickness / 2;
//        Interval intervalX = new Interval(-plugR, plugR + 1);
//        Interval intervalY = new Interval(-tabCutoutThickness / 2, tabCutoutThickness / 2);
//        Interval intervalZ = new Interval(-(1 + eps + tabCutoutDepth) / 2, (eps + tabCutoutDepth) / 2);
//        Brep cutoutBox = new Box(b, intervalX, intervalY, intervalZ).ToBrep();
//        var xf = Transform.Translation(0, ell, tabCutoutDepth / 2);
//        cutoutBox.Transform(xf);
//        return cutoutBox;
//    }
//    /**
//# Make a triangle for wedge profile with right angle at the origin
//# height(along zaxis) = plug_wedge_h, the height of the wedge
//# length(along yaxis) = plug_wedge_l+plug_tab_thickness, the lenght increases as the gap nears the center
//    p0=rg.Point3d(0,0,0)
//p1=rg.Point3d(0,0,-plug_wedge_h)
//p2=rg.Point3d(0,plug_wedge_l+plug_tab_thickness,0)
//extrudePt=rg.Point3d(plug_wedge_w,0,0) #point to extrude the profile out to

//#build the wedge profile in the 4th quad, extrude the profile to plug_wedge_w, and then re-center the width on the xaxis
//wedge=rs.MoveObject(
//    rs.ExtrudeCurveStraight(
//        rs.AddPolyline([p1, p0, p2, p1]),
//    p0,extrudePt),
//rg.Vector3d(-plug_wedge_w/2,0,0))
//    */

//    public Brep rawWedge()
//    {
//        //Make a triangle for wedge profile with right angle at the origin
//        //height(along zaxis) = plug_wedge_h, the height of the wedge
//        //length(along yaxis) = plug_wedge_l+plug_tab_thickness, the lenght increases as the gap nears the center
//        Point3d p0 = new Point3d(0, 0, 0);
//        Point3d p1 = new Point3d(0, 0, -wedgeHeight);
//        Point3d p2 = new Point3d(0, wedgeLength + tabThickness, 0);
//        Vector3d extrudePath = new Vector3d(wedgeWidth, 0, 0); //path to extrude along = width of wedge

//        //build the wedge profile in the 4th quad
//        PolylineCurve profile = new PolylineCurve(new Point3dList(p1, p0, p2, p1));
//        //extrude the profile to plug_wedge_w
//        Brep wedge = Extrusion.Create(profile, wedgeWidth, true).ToBrep();
//        //and then re-center the width on the xaxis
//        var xf = Transform.Translation(wedgeWidth / 2, 0, 0);
//        wedge.Transform(xf);
//        return wedge;
//    }
//}

namespace MyPlugin.UtilityComps
{
    public class PositivePlugComponent : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the PositivePlugComponent class.
        /// </summary>
        public PositivePlugComponent()
          : base("PositivePlug", "PP",
              "Description",
              "MyPlugin", "Utilites")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("wedgeHeight", "WH", "height of tab wedge", GH_ParamAccess.item, 4.0);
            pManager.AddNumberParameter("wedgeLength","WL", "length of tab wedge", GH_ParamAccess.item,2.0);
            pManager.AddNumberParameter("wedgeWidth","WW","width of tab wedge",GH_ParamAccess.item,2.0);
            pManager.AddNumberParameter("tabCutoutThickness","CT","space between tab and center of plug",GH_ParamAccess.item,1.40);
            pManager.AddNumberParameter("tabCutoutDepth","CD","how far down do the tabs go",GH_ParamAccess.item,10.0);
            pManager.AddNumberParameter("tabThickness","TT","thickness of both tabs",GH_ParamAccess.item,2.76);
            pManager.AddNumberParameter("taperfactor", "TF", "how much smaller the tapered radius is to the plug radius", GH_ParamAccess.item,0.70);
            pManager.AddNumberParameter("TaperLength","TL", "distance between the smaller radius and main body cyllinder", GH_ParamAccess.item,5.0);
            pManager.AddNumberParameter("socketCylDia", "SD","Diameter of socket",GH_ParamAccess.item,15.00);
            pManager.AddNumberParameter("socketWallThickness","ST", "thickness from the outer diameter of the socket inward", GH_ParamAccess.item,2.00);
            pManager.AddNumberParameter("connectionPlay","CP", "an adjustable amount to tweak the snuggness of the plugs fit into the socket", GH_ParamAccess.item,0.20);
            pManager.AddNumberParameter("socketLength","SL", "height of the socket",GH_ParamAccess.item,7.0);
            pManager.AddNumberParameter("bodyOverlap","BO", "how much the plug and socket overlap", GH_ParamAccess.item,7.0);
            pManager.AddNumberParameter("lengthOverage", "LO", "how much excess hangover between the plug and socket", GH_ParamAccess.item,7.0);
            pManager.AddNumberParameter("eps","E", "an adjusment value default 0.01", GH_ParamAccess.item,0.1);
            pManager.AddPlaneParameter("base", "B", "base xy plane to build upon", GH_ParamAccess.item, Plane.WorldXY);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGeometryParameter("plugBody", "G", "cylinder brep to be subtracted from the surface", GH_ParamAccess.item);
            pManager.AddGeometryParameter("PBOneWedge", "P1","first wedge", GH_ParamAccess.item);
            pManager.AddGeometryParameter("PBTwoWedge", "P2","second wedge", GH_ParamAccess.item);
            pManager.AddGeometryParameter("wedge1", "W1", "first cutout", GH_ParamAccess.item);
            pManager.AddGeometryParameter("wedge2", "W2", "seocnd cutout", GH_ParamAccess.item);
           // pManager.AddIntegerParameter("length", "L", "seocnd cutout", GH_ParamAccess.item);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double wedgeHeight = 0;
            double wedgeLength = 0;
            double wedgeWidth = 0;
            double tabCutoutThickness = 0;
            double tabCutoutDepth = 0;
            double tabThickness = 0;
            double taperfactor = 0;
            double plugTaperLength = 0;
            double socketCylDia = 0;
            double socketWallThickness = 0;
            double connectionPlay = 0;
            double socketLength = 0;
            double bodyOverlap = 0;
            double lengthOverage = 0;
            double eps = 0;
            Plane b = Plane.Unset;
            
            if(!DA.GetData(0, ref wedgeHeight)) return;
            if(!DA.GetData(1, ref wedgeLength)) return;
            if(!DA.GetData(2, ref wedgeWidth)) return;
            if(!DA.GetData(3, ref tabCutoutThickness)) return;
            if(!DA.GetData(4, ref tabCutoutDepth)) return;
            if(!DA.GetData(5, ref tabThickness)) return;
            if(!DA.GetData(6, ref taperfactor)) return;
            if(!DA.GetData(7, ref plugTaperLength)) return;
            if(!DA.GetData(8, ref socketCylDia)) return;
            if(!DA.GetData(9, ref socketWallThickness)) return;
            if(!DA.GetData(10, ref connectionPlay)) return;
            if(!DA.GetData(11, ref socketLength)) return;
            if(!DA.GetData(12, ref bodyOverlap)) return;
            if(!DA.GetData(13, ref lengthOverage)) return;
            if(!DA.GetData(14, ref eps)) return;
            if(!DA.GetData(15, ref b)) return;
            
            double plugR = (socketCylDia / 2) - socketWallThickness - connectionPlay; //radius of the plug shuld fit snug into the socket
            double taperedR = taperfactor * plugR;
            //main body height fits entirely into the socket plus some overlap and excess length, minus the taper length
            double bodyHeight = socketLength + bodyOverlap - plugTaperLength + lengthOverage;
            
            /**Send data to objects for creation*/
            PlugBody plugBody = new PlugBody(taperedR, plugTaperLength, plugR, bodyHeight); //used to create the plug body
            PlugTabs plugTabs = new PlugTabs(wedgeHeight, wedgeLength, wedgeWidth, tabCutoutThickness, tabCutoutDepth, tabThickness); //used to create cuoutBoxes and Wedges

            /** create the overall plug body*/
            Brep plugBodyBrep = plugBody.plugBodyGeo();
            DA.SetData(0, plugBodyBrep); //send the plug body geometry to "plugBody" output for debugging

            /**Create and place cutout boxes to create the tabs
             * The difference of the boxes will be taken from the plugbody
             * Cutout boxes should be equi distance from the center of the plug
             */
            // var mf = Transform.Mirror(Plane.WorldZX); //create mirror transformation vector matrix
            Brep[] cutouts= new Brep[2] { plugTabs.cutoutBox(plugR, b, eps), plugTabs.cutoutBox(plugR, b, eps) }; //create two new cutout boxes
            plugBodyBrep = cutoutBoxes(cutouts, plugBodyBrep);
        
            
            /**Create and place wedges onto outside of plug body
             * Vertical Face should be parrallel to cutout box long face*/
            Brep[] wedges = new Brep[2] { plugTabs.rawWedge(), plugTabs.rawWedge() }; //Create 2 wedges
           ;
            plugBodyBrep = addWedges(wedges, plugR, tabThickness, plugTaperLength, plugBodyBrep);
            DA.SetData(2, plugBodyBrep);
            
        }

        public Brep addWedges(Brep[] wedges,double plugR,double tabThickness, double plugTaperLength, Brep plugBodyBrep)
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
        public Brep cutoutBoxes(Brep[] cutouts, Brep plugBodyBrep)
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
                return Properties.Resources.positive_plug_icon;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("ED89E2C6-6DBE-4208-9B52-70B1A2B53814"); }
        }
    }
}

//using System;
//using System.Collections;
//using System.Collections.Generic;

//using Grasshopper.Kernel;
//using Rhino;
//using Rhino.Collections;
//using Rhino.Geometry;
//using Rhino.Geometry.Collections;
//using Rhino.PlugIns;
//using Rhino.Render.ChangeQueue;
//using Rhino.Render.PostEffects;

//public class PlugBody
//{
//    public double taperedR { get; set; }
//    public double plugTaperLength { get; set; }
//    public double plugR { get; set; }
//    public double bodyHeight { get; set; }
//    public PlugBody(double taperedR, double plugTaperLength, double plugR, double bodyHeight)
//    {
//        this.taperedR = taperedR;
//        this.plugTaperLength = plugTaperLength;
//        this.plugR = plugR;
//        this.bodyHeight = bodyHeight;
//    }
//    public Brep plugBodyGeo()
//    {
//        Point3d p0 = new Point3d(0, 0, 0);
//        Point3d p1 = new Point3d(taperedR, 0, 0); //origin to taper radius
//        Point3d p2 = new Point3d(plugR, 0, plugTaperLength); //transition from taper to body
//        Point3d p3 = new Point3d(plugR, 0, plugTaperLength + bodyHeight); //to top of plug
//        Point3d p4 = new Point3d(0, 0, plugTaperLength + bodyHeight); //top on axis

//        Polyline polyline = new Polyline(new Point3dList(p0, p1, p2, p3));
//        Line axis = new Line(p0, p4);
//        Brep geo = RevSurface.Create(polyline, axis).ToBrep();
//        return geo.CapPlanarHoles(0.001);
//    }
//}

//public class PlugTabs
//{
//    public double wedgeHeight { get; set; }
//    public double wedgeLength { get; set; }
//    public double wedgeWidth { get; set; }
//    public double tabCutoutThickness { get; set; }
//    public double tabCutoutDepth { get; set; }
//    public double tabThickness { get; set; }

//    public PlugTabs(double wedgeHeight, double wedgeLength, double wedgeWidth, double tabCutoutThickness, double tabCutoutDepth, double tabThickness)
//    {
//        this.wedgeHeight = wedgeHeight;
//        this.wedgeLength = wedgeLength;
//        this.wedgeWidth = wedgeWidth;
//        this.tabCutoutThickness = tabCutoutThickness;
//        this.tabCutoutDepth = tabCutoutDepth;
//        this.tabThickness = tabThickness;
//    }
//    public Brep cutoutBox(double plugR, Plane b, double eps)
//    {
//        double ell = plugR - tabThickness - tabCutoutThickness / 2;
//        Interval intervalX = new Interval(-plugR, plugR + 1);
//        Interval intervalY = new Interval(-tabCutoutThickness / 2, tabCutoutThickness / 2);
//        Interval intervalZ = new Interval(-(1 + eps + tabCutoutDepth) / 2, (eps + tabCutoutDepth) / 2);
//        Brep cutoutBox = new Box(b, intervalX, intervalY, intervalZ).ToBrep();
//        var xf = Transform.Translation(0, ell, tabCutoutDepth / 2);
//        cutoutBox.Transform(xf);
//        return cutoutBox;
//    }
//    /**
//# Make a triangle for wedge profile with right angle at the origin
//# height(along zaxis) = plug_wedge_h, the height of the wedge
//# length(along yaxis) = plug_wedge_l+plug_tab_thickness, the lenght increases as the gap nears the center
//    p0=rg.Point3d(0,0,0)
//p1=rg.Point3d(0,0,-plug_wedge_h)
//p2=rg.Point3d(0,plug_wedge_l+plug_tab_thickness,0)
//extrudePt=rg.Point3d(plug_wedge_w,0,0) #point to extrude the profile out to

//#build the wedge profile in the 4th quad, extrude the profile to plug_wedge_w, and then re-center the width on the xaxis
//wedge=rs.MoveObject(
//    rs.ExtrudeCurveStraight(
//        rs.AddPolyline([p1, p0, p2, p1]),
//    p0,extrudePt),
//rg.Vector3d(-plug_wedge_w/2,0,0))
//    */

//    public Brep rawWedge()
//    {
//        //Make a triangle for wedge profile with right angle at the origin
//        //height(along zaxis) = plug_wedge_h, the height of the wedge
//        //length(along yaxis) = plug_wedge_l+plug_tab_thickness, the lenght increases as the gap nears the center
//        Point3d p0 = new Point3d(0, 0, 0);
//        Point3d p1 = new Point3d(0, 0, -wedgeHeight);
//        Point3d p2 = new Point3d(0, wedgeLength + tabThickness, 0);
//        Vector3d extrudePath = new Vector3d(wedgeWidth, 0, 0); //path to extrude along = width of wedge

//        //build the wedge profile in the 4th quad
//        PolylineCurve profile = new PolylineCurve(new Point3dList(p1, p0, p2, p1));
//        //extrude the profile to plug_wedge_w
//        Brep wedge = Extrusion.Create(profile, wedgeWidth, true).ToBrep();
//        //and then re-center the width on the xaxis
//        var xf = Transform.Translation(wedgeWidth / 2, 0, 0);
//        wedge.Transform(xf);
//        return wedge;
//    }
//}

//namespace MyPlugin.UtilityComps
//{
//    public class PositivePlugComponent : GH_Component
//    {
//        /// <summary>
//        /// Initializes a new instance of the PositivePlugComponent class.
//        /// </summary>
//        public PositivePlugComponent()
//          : base("PositivePlug", "PP",
//              "Description",
//              "MyPlugin", "Utilites")
//        {
//        }

//        /// <summary>
//        /// Registers all the input parameters for this component.
//        /// </summary>
//        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
//        {
//            pManager.AddNumberParameter("wedgeHeight", "WH", "height of tab wedge", GH_ParamAccess.item, 4.0);
//            pManager.AddNumberParameter("wedgeLength", "WL", "length of tab wedge", GH_ParamAccess.item, 2.0);
//            pManager.AddNumberParameter("wedgeWidth", "WW", "width of tab wedge", GH_ParamAccess.item, 2.0);
//            pManager.AddNumberParameter("tabCutoutThickness", "CT", "space between tab and center of plug", GH_ParamAccess.item, 1.40);
//            pManager.AddNumberParameter("tabCutoutDepth", "CD", "how far down do the tabs go", GH_ParamAccess.item, 10.0);
//            pManager.AddNumberParameter("tabThickness", "TT", "thickness of both tabs", GH_ParamAccess.item, 2.76);
//            pManager.AddNumberParameter("taperfactor", "TF", "how much smaller the tapered radius is to the plug radius", GH_ParamAccess.item, 0.70);
//            pManager.AddNumberParameter("TaperLength", "TL", "distance between the smaller radius and main body cyllinder", GH_ParamAccess.item, 5.0);
//            pManager.AddNumberParameter("socketCylDia", "SD", "Diameter of socket", GH_ParamAccess.item, 15.00);
//            pManager.AddNumberParameter("socketWallThickness", "ST", "thickness from the outer diameter of the socket inward", GH_ParamAccess.item, 2.00);
//            pManager.AddNumberParameter("connectionPlay", "CP", "an adjustable amount to tweak the snuggness of the plugs fit into the socket", GH_ParamAccess.item, 0.20);
//            pManager.AddNumberParameter("socketLength", "SL", "height of the socket", GH_ParamAccess.item, 7.0);
//            pManager.AddNumberParameter("bodyOverlap", "BO", "how much the plug and socket overlap", GH_ParamAccess.item, 7.0);
//            pManager.AddNumberParameter("lengthOverage", "LO", "how much excess hangover between the plug and socket", GH_ParamAccess.item, 7.0);
//            pManager.AddNumberParameter("eps", "E", "an adjusment value default 0.01", GH_ParamAccess.item, 0.1);
//            pManager.AddPlaneParameter("base", "B", "base xy plane to build upon", GH_ParamAccess.item, Plane.WorldXY);

//        }

//        /// <summary>
//        /// Registers all the output parameters for this component.
//        /// </summary>
//        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
//        {
//            pManager.AddGeometryParameter("plugBody", "G", "cylinder brep to be subtracted from the surface", GH_ParamAccess.item);
//            pManager.AddGeometryParameter("PBOneWedge", "P1", "first wedge", GH_ParamAccess.item);
//            pManager.AddGeometryParameter("PBTwoWedge", "P2", "second wedge", GH_ParamAccess.item);
//            pManager.AddGeometryParameter("wedge1", "W1", "first cutout", GH_ParamAccess.item);
//            pManager.AddGeometryParameter("wedge2", "W2", "seocnd cutout", GH_ParamAccess.item);
//            // pManager.AddIntegerParameter("length", "L", "seocnd cutout", GH_ParamAccess.item);

//        }

//        /// <summary>
//        /// This is the method that actually does the work.
//        /// </summary>
//        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
//        protected override void SolveInstance(IGH_DataAccess DA)
//        {
//            double wedgeHeight = 0;
//            double wedgeLength = 0;
//            double wedgeWidth = 0;
//            double tabCutoutThickness = 0;
//            double tabCutoutDepth = 0;
//            double tabThickness = 0;
//            double taperfactor = 0;
//            double plugTaperLength = 0;
//            double socketCylDia = 0;
//            double socketWallThickness = 0;
//            double connectionPlay = 0;
//            double socketLength = 0;
//            double bodyOverlap = 0;
//            double lengthOverage = 0;
//            double eps = 0;
//            Plane b = Plane.Unset;

//            DA.GetData(0, ref wedgeHeight);
//            DA.GetData(1, ref wedgeLength);
//            DA.GetData(2, ref wedgeWidth);
//            DA.GetData(3, ref tabCutoutThickness);
//            DA.GetData(4, ref tabCutoutDepth);
//            DA.GetData(5, ref tabThickness);
//            DA.GetData(6, ref taperfactor);
//            DA.GetData(7, ref plugTaperLength);
//            DA.GetData(8, ref socketCylDia);
//            DA.GetData(9, ref socketWallThickness);
//            DA.GetData(10, ref connectionPlay);
//            DA.GetData(11, ref socketLength);
//            DA.GetData(12, ref bodyOverlap);
//            DA.GetData(13, ref lengthOverage);
//            DA.GetData(14, ref eps);
//            DA.GetData(15, ref b);

//            double plugR = (socketCylDia / 2) - socketWallThickness - connectionPlay; //radius of the plug shuld fit snug into the socket
//            double taperedR = taperfactor * plugR;
//            //main body height fits entirely into the socket plus some overlap and excess length, minus the taper length
//            double bodyHeight = socketLength + bodyOverlap - plugTaperLength + lengthOverage;

//            /**Send data to objects for creation*/
//            PlugBody plugBody = new PlugBody(taperedR, plugTaperLength, plugR, bodyHeight); //used to create the plug body
//            PlugTabs plugTabs = new PlugTabs(wedgeHeight, wedgeLength, wedgeWidth, tabCutoutThickness, tabCutoutDepth, tabThickness); //used to create cuoutBoxes and Wedges

//            /** create the overall plug body*/
//            Brep plugBodyBrep = plugBody.plugBodyGeo();
//            DA.SetData(0, plugBodyBrep); //send the plug body geometry to "plugBody" output for debugging

//            /**Create and place cutout boxes to create the tabs
//             * The difference of the boxes will be taken from the plugbody
//             * Cutout boxes should be equi distance from the center of the plug
//             */
//            // var mf = Transform.Mirror(Plane.WorldZX); //create mirror transformation vector matrix
//            Brep[] cutouts = new Brep[2] { plugTabs.cutoutBox(plugR, b, eps), plugTabs.cutoutBox(plugR, b, eps) }; //create two new cutout boxes
//            plugBodyBrep = cutoutBoxes(cutouts, plugBodyBrep);
//            //DA.SetData(3, cutouts[0]); //send the cutout geometry to outparamter "cutout1" for debugging
//            //cutouts[1].Transform(mf); //mirror the second cutout to the other side of the plug body
//            //    //DA.SetData(4, cutouts[1]); //send the cutout geometry to outparamter "cutout2" for debugging

//            ///**Cutout the cutout boxes from the plug body to create tabs*/
//            //plugBodyBrep = Brep.CreateBooleanDifference(plugBodyBrep, cutouts[0], 0.01)[0]; //plugBody-cutout[0]
//            //    //DA.SetData(1, plugBodyBrep); //Send plugBody with one cutout to "PBOneCut" for debugging 
//            /////Because of the mirror, the polyline orientation of cutout[1] can become flipped from counter-clockwise to clockwise
//            /////Boolean operations can become flipped when their orientation is flipped
//            /////Read more: https://discourse.mcneel.com/t/rhino-geometry-brep-createbooleandifference-not-creating-difference/116008
//            /////Solution: https://discourse.mcneel.com/t/rhino-geometry-brep-createbooleandifference-not-creating-difference/116008
//            /////Without these 3 lines, the difference results in cutout-plugBody, not plugbody-cutout
//            //cutouts[1].Faces.SplitKinkyFaces(RhinoMath.DefaultAngleTolerance, true);
//            //if (BrepSolidOrientation.Inward == cutouts[1].SolidOrientation)
//            //{
//            //    cutouts[1].Flip();
//            //}
//            //plugBodyBrep = Brep.CreateBooleanDifference(plugBodyBrep, cutouts[1], 0.01)[0]; //Cutout the second cutout box. Need Split KinkyFaces and Flip Orientation Lines! 
//            //    //DA.SetData(2, plugBodyBrep); //send the plug body with both cutouts to "PBTwoCut" for debugging


//            /**Create and place wedges onto outside of plug body
//             * Vertical Face should be parrallel to cutout box long face*/
//            Brep[] wedges = new Brep[2] { plugTabs.rawWedge(), plugTabs.rawWedge() }; //Create 2 wedges
//            //var xf = Transform.Translation(0, plugR - tabThickness, plugTaperLength); //create transformation matrix for both wedges
//            //wedges[0].Transform(xf); //move wedges into place
//            //DA.SetData(3, wedges[0]);
//            //wedges[1].Transform(xf);
//            //wedges[1].Transform(Transform.Mirror(Plane.WorldZX)); //Flip one wedge to the other side of the plug
//            //DA.SetData(4, wedges[1]);
//            ////Brep[] thingsToUnion = new Brep[2] { plugBodyBrep, wedges[0] };
//            ////plugBodyBrep = Brep.CreateBooleanUnion(thingsToUnion, 0.01)[0];
//            //DA.SetData(1,plugBodyBrep);
//            //wedges[1].Faces.SplitKinkyFaces(RhinoMath.DefaultAngleTolerance, true);
//            //if (BrepSolidOrientation.Inward == wedges[1].SolidOrientation)
//            //{
//            //    wedges[1].Flip();
//            //}
//            //Brep[] thingsToUnion = new Brep[3] { plugBodyBrep, wedges[0], wedges[1] };
//            //thingsToUnion[1] = wedges[1];
//            //plugBodyBrep = Brep.CreateBooleanUnion(addWedges, 0.01)[0];
//            plugBodyBrep = addWedges(wedges, plugR, tabThickness, plugTaperLength, plugBodyBrep);
//            DA.SetData(2, plugBodyBrep);

//        }

//        public Brep addWedges(Brep[] wedges, double plugR, double tabThickness, double plugTaperLength, Brep plugBodyBrep)
//        {
//            /**Create and place wedges onto outside of plug body
//             * Vertical Face should be parrallel to cutout box long face*/
//            var xf = Transform.Translation(0, plugR - tabThickness, plugTaperLength); //create transformation matrix for both wedges
//            wedges[0].Transform(xf); //move wedges into place

//            wedges[1].Transform(xf);
//            wedges[1].Transform(Transform.Mirror(Plane.WorldZX)); //Flip one wedge to the other side of the plug

//            //Brep[] thingsToUnion = new Brep[2] { plugBodyBrep, wedges[0] };
//            //plugBodyBrep = Brep.CreateBooleanUnion(thingsToUnion, 0.01)[0];

//            wedges[1].Faces.SplitKinkyFaces(RhinoMath.DefaultAngleTolerance, true);
//            if (BrepSolidOrientation.Inward == wedges[1].SolidOrientation)
//            {
//                wedges[1].Flip();
//            }
//            Brep[] thingsToUnion = new Brep[3] { plugBodyBrep, wedges[0], wedges[1] };

//            //thingsToUnion[1] = wedges[1];
//            plugBodyBrep = Brep.CreateBooleanUnion(thingsToUnion, 0.01)[0];
//            return plugBodyBrep;
//        }
//        public Brep cutoutBoxes(Brep[] cutouts, Brep plugBodyBrep)
//        {
//            /**Create and place cutout boxes to create the tabs
//                 * The difference of the boxes will be taken from the plugbody
//                 * Cutout boxes should be equi distance from the center of the plug
//                 */
//            var mf = Transform.Mirror(Plane.WorldZX); //create mirror transformation vector matrix
//                                                      //Brep[] cutouts = new Brep[2] { plugTabs.cutoutBox(plugR, b, eps), plugTabs.cutoutBox(plugR, b, eps) }; //create two new cutout boxes
//                                                      //DA.SetData(3, cutouts[0]); //send the cutout geometry to outparamter "cutout1" for debugging
//            cutouts[1].Transform(mf); //mirror the second cutout to the other side of the plug body
//                                      //DA.SetData(4, cutouts[1]); //send the cutout geometry to outparamter "cutout2" for debugging

//            /**Cutout the cutout boxes from the plug body to create tabs*/
//            plugBodyBrep = Brep.CreateBooleanDifference(plugBodyBrep, cutouts[0], 0.01)[0]; //plugBody-cutout[0]
//                                                                                            //DA.SetData(1, plugBodyBrep); //Send plugBody with one cutout to "PBOneCut" for debugging 
//            ///Because of the mirror, the polyline orientation of cutout[1] can become flipped from counter-clockwise to clockwise
//            ///Boolean operations can become flipped when their orientation is flipped
//            ///Read more: https://discourse.mcneel.com/t/rhino-geometry-brep-createbooleandifference-not-creating-difference/116008
//            ///Solution: https://discourse.mcneel.com/t/rhino-geometry-brep-createbooleandifference-not-creating-difference/116008
//            ///Without these 3 lines, the difference results in cutout-plugBody, not plugbody-cutout
//            cutouts[1].Faces.SplitKinkyFaces(RhinoMath.DefaultAngleTolerance, true);
//            if (BrepSolidOrientation.Inward == cutouts[1].SolidOrientation)
//            {
//                cutouts[1].Flip();
//            }
//            plugBodyBrep = Brep.CreateBooleanDifference(plugBodyBrep, cutouts[1], 0.01)[0]; //Cutout the second cutout box. Need Split KinkyFaces and Flip Orientation Lines! 
//            return plugBodyBrep; //DA.SetData(2, plugBodyBrep); //send the plug body with both cutouts to "PBTwoCut" for debugging

//        }

//        /// <summary>
//        /// Provides an Icon for the component.
//        /// </summary>
//        protected override System.Drawing.Bitmap Icon
//        {
//            get
//            {
//                //You can add image files to your project resources and access them like this:
//                // return Resources.IconForThisComponent;
//                return null;
//            }
//        }

//        /// <summary>
//        /// Gets the unique ID for this component. Do not change this ID after release.
//        /// </summary>
//        public override Guid ComponentGuid
//        {
//            get { return new Guid("ED89E2C6-6DBE-4208-9B52-70B1A2B53814"); }
//        }
//    }
//}