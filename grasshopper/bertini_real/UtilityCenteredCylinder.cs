using System;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace bertini_real
{
    /// <summary>
    /// A closed (capped) cylinder/prism centered on a plane's origin and running along its Z axis
    /// from -Length/2 to +Length/2.  Sides picks the cross-section: 1 = a true round cylinder,
    /// 3 = triangular prism, 4 = square, 5 = pentagon, ...  (Sides = 2 is degenerate.)  Radius is
    /// the circumradius (distance to the polygon vertices).  Handy as a connector / boolean cutter
    /// -- centered so it straddles its placement point, closed so mesh booleans bite.
    /// </summary>
    public class UtilityCenteredCylinder : GH_Component
    {
        public UtilityCenteredCylinder()
          : base("Centered Closed Cylinder", "CenCyl",
                 "A capped cylinder/prism centered on the plane origin along its Z axis. Sides: 1 = round, 3+ = N-gon prism.",
                 "bertini_real", "Utility")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Radius", "R", "Radius (circumradius for polygons -- distance to the vertices)", GH_ParamAccess.item, 5.0);
            pManager.AddNumberParameter("Length", "L", "Length (centered on the plane origin)", GH_ParamAccess.item, 100.0);
            pManager.AddIntegerParameter("Sides", "N", "Cross-section: 1 = round cylinder, 3 = triangle, 4 = square, 5 = pentagon, ...", GH_ParamAccess.item, 1);
            pManager.AddPlaneParameter("Plane", "P", "Center plane; the cylinder runs along its Z axis", GH_ParamAccess.item, Plane.WorldXY);
            Params.Input[0].Optional = true;
            Params.Input[1].Optional = true;
            Params.Input[2].Optional = true;
            Params.Input[3].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddBrepParameter("Cylinder", "C", "Closed cylinder/prism Brep", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double radius = 5.0;
            double length = 100.0;
            int sides = 1;
            Plane plane = Plane.WorldXY;
            DA.GetData(0, ref radius);
            DA.GetData(1, ref length);
            DA.GetData(2, ref sides);
            DA.GetData(3, ref plane);

            if (radius <= 0 || length <= 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Radius and Length must be positive.");
                return;
            }
            if (sides < 1 || sides == 2)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Sides must be 1 (round) or >= 3 (polygon); 2 is degenerate.");
                return;
            }

            // base plane is half the length below the center, so the result is centered
            var bottom = new Plane(plane.Origin - plane.ZAxis * (length / 2.0), plane.XAxis, plane.YAxis);

            Brep brep;
            if (sides == 1)
            {
                var cylinder = new Cylinder(new Circle(bottom, radius), length);
                brep = cylinder.ToBrep(true, true);   // cap both ends -> closed solid
            }
            else
            {
                // regular n-gon profile in the base plane, extruded by Length and capped
                var pts = new Point3d[sides + 1];
                for (int i = 0; i < sides; i++)
                {
                    double a = 2.0 * Math.PI * i / sides;
                    pts[i] = bottom.PointAt(radius * Math.Cos(a), radius * Math.Sin(a));
                }
                pts[sides] = pts[0];

                var profile = new Polyline(pts).ToPolylineCurve();
                Extrusion extrusion = Extrusion.Create(profile, length, true);
                if (extrusion == null)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Failed to build the prism.");
                    return;
                }
                brep = extrusion.ToBrep();
            }

            DA.SetData(0, brep);
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("lego.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("8A3F9C21-6D74-4E08-B5A2-1C7E0F934D6B"); }
        }
    }
}
