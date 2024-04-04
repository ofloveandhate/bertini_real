using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace MyPlugin.UtilityComps
{
    public class NegativePlugComponent : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the NegativePlugComponent class.
        /// </summary>
        public NegativePlugComponent()
          : base("NegativePlug", "NegPlug",
              "Creates hole to remove from the piece and positive plug to allow threading a wire",
              "MyPlugin", "Utilites")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Wire Hole Diameter", "D", "diameter of hole required for a wire", GH_ParamAccess.item, 7.0);
            pManager.AddNumberParameter("Plug Length Factor", "PF", "multiplier for length of plug neg", GH_ParamAccess.item, 3.40);
            pManager.AddNumberParameter("Socket Length", "SL", "height of the socket", GH_ParamAccess.item, 7.0);
            pManager.AddNumberParameter("Length Overage", "LO", "how much excess hangover between the plug and socket", GH_ParamAccess.item, 7.0);
            pManager.AddNumberParameter("Body Overlap", "BO", "how much the plug and socket overlap", GH_ParamAccess.item, 7.0);
            pManager.AddNumberParameter("eps", "E", "an adjusment value default 0.01", GH_ParamAccess.item, 0.01);
            pManager.AddPlaneParameter("base", "B", "base xy plane to build upon", GH_ParamAccess.item, Plane.WorldXY);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGeometryParameter("geo", "G", "brep cylinder to be removed from the surface", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            ///define place holder vars
            ///These names can vary and are used within the code , but if you cnage their name ensure you change it throughout the code
            double wireHoleDia = 0;
            double plugFactor = 0;
            double socketLength = 0;
            double lengthOverage = 0;
            double bodyOverlap = 0;
            double eps = 0;
            Plane b = Plane.Unset;

            //retrieve inputs
            if (!DA.GetData(0, ref wireHoleDia)) return;
            if (!DA.GetData(1, ref plugFactor)) return;
            if (!DA.GetData(22, ref socketLength)) return;
            if (!DA.GetData(3, ref lengthOverage)) return;
            if (!DA.GetData(4, ref bodyOverlap))return;
            if(!DA.GetData(5, ref eps))return;
            if(!DA.GetData(6, ref b))return;

            /* Create a cylinder  */
            double cuttingR = wireHoleDia / 2; //set the radius
            double length = plugFactor * (socketLength + lengthOverage + bodyOverlap); //set the length
            
            ///Create the actual cyllinder in gh
            Brep cylinder = new Cylinder(new Circle(b, cuttingR), length).ToBrep(true, true); //create a cyllinder Brep with caps (can remove caps by changing true values) 
            var xf = Transform.Translation(0, 0, -eps - (length / 2)); //move it down to the center so transformations work easily later
            cylinder.Transform(xf);

            /*Send to output*/
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
                return Properties.Resources.negative_plug_icon;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("D6184726-4799-45BB-A812-7CD12FF66BCB"); }
        }
    }
}

/**The python code used
 * 
   """Creates hole to remove from the surface to insert plug
    Inputs:
        wire_hole_dia: diameter of hole required for a wire
        plug_neg_factor: multiplier for length of plug neg
        socket_length: lenght of socket
        eps: an adjusment amount
        base: base xy plane to build upon
        plug_length_overage: excess plug length
        plug_body_overlap: how much the body and surface overlap
    Output:
        geo: brep cylinder to be removed from the surface"""

    import rhinoscriptsyntax as rs
    import Rhino.Geometry as rg
    
    cutting_r = wire_hole_dia/2
    length = plug_neg_factor*(socket_length+plug_length_overage+plug_body_overlap)
    
    #create a cylinder, move down on z axis to center it on origin
    geo = rs.MoveObject(rs.AddCylinder(base,length,cutting_r),rg.Vector3d(0,0,-eps-length/2))

 */