using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;
using static Rhino.Render.TextureGraphInfo;

namespace MyPlugin.UtilityComps
{
/**
""creates the socket that is inserted into the surface
    Inputs:
        socket_cyl_dia: outer diameter of the socket
        socket_wall_thickness: thickness of socket from the socket_cyl_dia
        socket_length: length of socket
        eps: an adjusment value
        base: base xy plane to build upon
    Output:
        geo: socket brep"""

import rhinoscriptsyntax as rs
import Rhino.Geometry as rg

origin = rg.Point3d(0,0,0)
yaxis = rg.Vector3d(0,1,0)

#a cylinder centered at 0,0 with length in the -z direction
#create 2 cocentric cylinders, the inner cyllinder is offset from the outer cylinder by socket_wall_thickness
geo = rs.RotateObject(
    rs.BooleanDifference(
        rs.AddCylinder(base, socket_length,socket_cyl_dia/2),
        rs.MoveObject(
            rs.AddCylinder(base, socket_length+2*eps,socket_cyl_dia/2-socket_wall_thickness),
        rg.Vector3d(0,0,-eps))),
origin,180,yaxis)
*/
    public class PositiveSocketComponent : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public PositiveSocketComponent()
          : base("Positive Socket ", "PS",
              "Create the connector pieces for nodes.",
              "MyPlugin", "Utilites")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("SocketCylDia", "D", "Outer diameter of the socket", GH_ParamAccess.item, 15.00);
            pManager.AddNumberParameter("socket_wall_thickness", "T", "thickness of socket from the socket_cyl_dia", GH_ParamAccess.item, 2.00);
            pManager.AddNumberParameter("socket_length", "L", "length of positive socket default", GH_ParamAccess.item, 10.00);
            pManager.AddNumberParameter("eps", "E", "an adjusment value default 0.01", GH_ParamAccess.item, 0.01);
            pManager.AddPlaneParameter("base", "B", "base xy plane to build upon", GH_ParamAccess.item, Plane.WorldXY);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGeometryParameter("geo", "G", "cylinder brep to be subtracted from the surface", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            /**
            # a cylinder centered at 0,0 with length in the -z direction
            # create 2 cocentric cylinders, the inner cyllinder is offset from the outer cylinder by socket_wall_thickness
            geo = rs.RotateObject(
                rs.BooleanDifference(
                    rs.AddCylinder(base, socket_length, socket_cyl_dia / 2),
                    rs.MoveObject(
                        rs.AddCylinder(base, socket_length + 2 * eps, socket_cyl_dia / 2 - socket_wall_thickness),
                    rg.Vector3d(0, 0, -eps))),
            origin, 180, yaxis)
            */
            //define place holder vars
            double dia = 0;
            double socketThickness = 0;
            double socketLength = 0;
            double eps = 0;
            Plane b = Plane.Unset;

            //retrieve inputs allow to get data even if not set because they have defaults)
            if(!DA.GetData(0, ref dia))return;
            if(!DA.GetData(1, ref socketThickness))return;
            if(!DA.GetData(2, ref socketLength))return;
            if(!DA.GetData(3, ref eps))return;
            if(!DA.GetData(4, ref b))return;
            

            Brep innerCylinder = new Cylinder(new Circle(b, (dia / 2) - socketThickness), socketLength + 2 * eps).ToBrep(true,true);
            var xf = Transform.Translation(0, 0, -eps);
            Brep outerCylinder = new Cylinder(new Circle(b, dia / 2), socketLength).ToBrep(true, true);

            Brep[] socket = Brep.CreateBooleanDifference(outerCylinder, innerCylinder, 0.01);
            var rf = Transform.Rotation(Math.PI, Vector3d.YAxis, Point3d.Origin);
            socket[0].Transform(rf);
            DA.SetData(0, socket[0]);
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
                return Properties.Resources.positive_socket_icon;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("1F11D0EB-95A4-477D-8C4E-42772D29F55E"); }
        }
    }
}