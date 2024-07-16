using System;
using System.Collections.Generic;
using System.Drawing;
using System.Net.Sockets;
using System.Threading;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;
using static Rhino.Render.TextureGraphInfo;

/*<summary>
 * Build a Component for the Negative Socket Piece
 * This is an idividual component. This should be differenced from the piece
 * </summary>*/
namespace MyPlugin.UtilityComps
{
    public class SocketNegativeComponent : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public SocketNegativeComponent()
          : base("Negative Socket", "NegSoc",
              "Create a Negative Socket",
              "MyPlugin", "Utilites")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        /// <param name="pManager"></param> The worker which grabs all the parameters on the side of a component
        /// Every component listed in this class will appear on the left side of the component as a possible input
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Socket Diameter", "SD", "Outer diameter of the socket", GH_ParamAccess.item, 15.00);
            pManager.AddNumberParameter("Socket Factor","SF", "multiplier for negative socket length default", GH_ParamAccess.item,3.00);
            pManager.AddNumberParameter("Socket Wall Thickness","ST", "thickness of socket from the socket_cyl_dia default", GH_ParamAccess.item,2.00);
            pManager.AddNumberParameter("Socket Length","SL", "length of positive socket default", GH_ParamAccess.item,10.00);
            pManager.AddNumberParameter("eps","E","an adjusment value default 0.01",GH_ParamAccess.item,0.01);
            pManager.AddPlaneParameter("base", "B", "base xy plane to build upon", GH_ParamAccess.item,Plane.WorldXY);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            //Send out the negative socket geometry 
            pManager.AddGeometryParameter("Negative Socket", "G", "cylinder brep to be subtracted from the surface",GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
          
            //define place holder vars
            double socketDiameter = 0;
            double socketNegFactor = 0;
            double socketThickness = 0;
            double socketLength = 0;
            double eps = 0;
            Plane b = Plane.Unset;

            //retrieve inputs allow to get data even if not set because they have defaults)
            if(!DA.GetData(0, ref socketDiameter))return;
            if(!DA.GetData(1, ref socketNegFactor))return;
            if(!DA.GetData(2, ref socketThickness))return;
            if(!DA.GetData(3, ref socketLength))return;
            if(!DA.GetData(4, ref eps))return;
            if(!DA.GetData(5, ref b))return;

            //Point3d origin = rg.Point3d(0, 0, 0)
            //yaxis = rg.Vector3d(0, 1, 0)
            double length = socketNegFactor * socketLength+ eps;
            // length of cyllinder is the length of the socket plus some excess, scaled up the socket_neg_factor
            //create the the socket as a cyllinder from a base circle that has a radius of the inner socket 
            Brep cylinder = new Cylinder(new Circle(b, (socketDiameter / 2) - socketThickness), length).ToBrep(true,true);
            //rotate and move the cyllinder into correct position
            var xf = Transform.Translation(0, 0, -length/2);
            cylinder.Transform(xf);
            var rf = Transform.Rotation(Math.PI, Vector3d.YAxis, Point3d.Origin);
            cylinder.Transform(rf);
            DA.SetData(0, cylinder);

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
                return Properties.Resources.negative_socket_icon;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("DB3EEAFB-0064-4EB5-8FBA-A59604542026"); }
        }
    }
}

/**
 * Python Code
creates hole for socket to be place on the surface.
    Inputs:
socket_cyl_dia: outer diameter socket
socket_neg_factor: multiplier for negative socket length
        socket_wall_thickness: thickness of socket from the socket_cyl_dia
        socket_length: length of positive socket
        eps: an adjusment value
        base: base xy plane to build upon
    Output:
        geo: a cylinder brep to be subtracted from the surface"""
import rhinoscriptsyntax as rs
import Rhino.Geometry as rg

origin = rg.Point3d(0,0,0)
yaxis = rg.Vector3d(0,1,0)
#length of cyllinder is the length of the socket plus some excess, scaled up the socket_neg_factor
length = socket_neg_factor*socket_length+eps

#cylinder centered at origin
geo = rs.RotateObject(rs.MoveObject(rs.AddCylinder(base, length, socket_cyl_dia/2-socket_wall_thickness),rg.Vector3d(0,0,-length/2)),origin,180,yaxis)
*/