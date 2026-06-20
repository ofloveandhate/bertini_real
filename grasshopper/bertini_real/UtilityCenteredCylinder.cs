using System;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace bertini_real
{
    /// <summary>
    /// A closed (capped) cylinder centered on a plane's origin and running along its Z axis from
    /// -Length/2 to +Length/2.  Handy as a connector / boolean cutter: centered so it straddles
    /// the point it's placed at, and closed so mesh booleans actually bite -- without wiring a
    /// circle + cylinder + cap every time.
    /// </summary>
    public class UtilityCenteredCylinder : GH_Component
    {
        public UtilityCenteredCylinder()
          : base("Centered Closed Cylinder", "CenCyl",
                 "A capped cylinder centered on the plane origin, along its Z axis (-Length/2 .. +Length/2)",
                 "bertini_real", "Utility")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Radius", "R", "Cylinder radius", GH_ParamAccess.item, 5.0);
            pManager.AddNumberParameter("Length", "L", "Cylinder length (centered on the plane origin)", GH_ParamAccess.item, 100.0);
            pManager.AddPlaneParameter("Plane", "P", "Center plane; the cylinder runs along its Z axis", GH_ParamAccess.item, Plane.WorldXY);
            Params.Input[0].Optional = true;
            Params.Input[1].Optional = true;
            Params.Input[2].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddBrepParameter("Cylinder", "C", "Closed cylinder Brep", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double radius = 5.0;
            double length = 100.0;
            Plane plane = Plane.WorldXY;
            DA.GetData(0, ref radius);
            DA.GetData(1, ref length);
            DA.GetData(2, ref plane);

            if (radius <= 0 || length <= 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Radius and Length must be positive.");
                return;
            }

            // base circle is half the length below the center, so the cylinder is centered
            var bottom = new Plane(plane.Origin - plane.ZAxis * (length / 2.0), plane.XAxis, plane.YAxis);
            var cylinder = new Cylinder(new Circle(bottom, radius), length);
            Brep brep = cylinder.ToBrep(true, true);   // cap both ends -> closed solid

            DA.SetData(0, brep);
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("lego.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("8A3F9C21-6D74-4E08-B5A2-1C7E0F934D6B"); }
        }
    }
}
