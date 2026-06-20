using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper;
using Grasshopper.GUI.Gradient;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Expressions;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

namespace bertini_real
{
    /// <summary>
    /// Colors piece meshes by a scalar function of position, evaluated at every mesh vertex.  The
    /// Function is an expression in x, y, z (e.g. "x^2 + y^2 + z^2", "Sin(x)*z"); its values are
    /// normalized over all the meshes (or to an explicit Domain) and mapped through a color
    /// gradient onto the mesh vertex colors.
    ///
    /// Rhino meshes carry vertex colors, and Spread Pieces preserves them, so this can go before
    /// or after spreading; preview the output meshes to see the coloring.
    /// </summary>
    public class SurfaceColorByFunction : GH_Component
    {
        public SurfaceColorByFunction()
          : base("Color By Function", "ColorFn",
                 "Color piece meshes by a scalar function of x,y,z evaluated at each vertex",
                 "bertini_real", "Surface")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Meshes", "M", "Per-piece meshes to color (one branch per piece)", GH_ParamAccess.tree);
            pManager.AddTextParameter("Function", "F", "Scalar expression in x, y, z (e.g. x^2+y^2+z^2)", GH_ParamAccess.item, "z");
            pManager.AddColourParameter("Colours", "Cs", "Gradient stops (>= 2); default is a blue->red spectrum", GH_ParamAccess.list);
            pManager.AddIntervalParameter("Domain", "D", "Value range mapped onto the gradient; default = the data's min..max", GH_ParamAccess.item);
            Params.Input[1].Optional = true;
            Params.Input[2].Optional = true;
            Params.Input[3].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Meshes", "M", "Colored meshes (vertex colors set), same tree structure", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Values", "V", "Per-vertex function values, parallel to each mesh's vertices", GH_ParamAccess.tree);
            pManager.AddIntervalParameter("Domain", "D", "The value range used for the gradient", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<GH_Mesh> meshes;
            if (!DA.GetDataTree(0, out meshes)) return;

            string function = "z";
            DA.GetData(1, ref function);
            if (string.IsNullOrWhiteSpace(function))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Function expression is empty.");
                return;
            }

            var colours = new List<Color>();
            DA.GetDataList(2, colours);

            Interval domain = Interval.Unset;
            bool hasDomain = DA.GetData(3, ref domain);

            // pass 1: evaluate the function at every vertex, tracking the global min/max
            var entries = new List<(GH_Path path, Mesh mesh, double[] vals)>();
            double gmin = double.MaxValue, gmax = double.MinValue;
            var parser = new GH_ExpressionParser();

            for (int b = 0; b < meshes.PathCount; b++)
            {
                GH_Path path = meshes.get_Path(b);
                foreach (var goo in meshes.get_Branch(path))
                {
                    if (!(goo is GH_Mesh gm) || gm.Value == null) continue;
                    Mesh mesh = gm.Value;
                    int n = mesh.Vertices.Count;
                    var vals = new double[n];
                    for (int i = 0; i < n; i++)
                    {
                        Point3d p = mesh.Vertices[i];
                        parser.AddVariable("x", p.X);
                        parser.AddVariable("y", p.Y);
                        parser.AddVariable("z", p.Z);
                        double v;
                        try { v = parser.Evaluate(function)._Double; }
                        catch (Exception e)
                        {
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Could not evaluate function: " + e.Message);
                            return;
                        }
                        vals[i] = v;
                        if (v < gmin) gmin = v;
                        if (v > gmax) gmax = v;
                    }
                    entries.Add((path, mesh, vals));
                }
            }

            if (entries.Count == 0) return;

            double lo = hasDomain ? domain.Min : gmin;
            double hi = hasDomain ? domain.Max : gmax;
            double span = hi - lo;
            if (Math.Abs(span) < 1e-15) span = 1.0;

            GH_Gradient gradient = (colours != null && colours.Count >= 2) ? BuildGradient(colours) : DefaultGradient();

            // pass 2: color each mesh and emit
            var outMesh = new DataTree<Mesh>();
            var outValues = new DataTree<double>();

            foreach (var (path, mesh, vals) in entries)
            {
                Mesh dup = mesh.DuplicateMesh();
                dup.VertexColors.Clear();
                for (int i = 0; i < vals.Length; i++)
                {
                    double t = (vals[i] - lo) / span;
                    if (t < 0) t = 0; else if (t > 1) t = 1;
                    dup.VertexColors.Add(gradient.ColourAt(t));
                }
                outMesh.Add(dup, path);
                outValues.AddRange(vals, path);
            }

            DA.SetDataTree(0, outMesh);
            DA.SetDataTree(1, outValues);
            DA.SetData(2, new Interval(lo, hi));
        }

        private static GH_Gradient BuildGradient(List<Color> colours)
        {
            var g = new GH_Gradient();
            int n = colours.Count;
            for (int i = 0; i < n; i++)
                g.AddGrip((double)i / (n - 1), colours[i]);
            return g;
        }

        private static GH_Gradient DefaultGradient()
        {
            var g = new GH_Gradient();
            g.AddGrip(0.00, Color.FromArgb(0, 0, 200));
            g.AddGrip(0.25, Color.FromArgb(0, 200, 200));
            g.AddGrip(0.50, Color.FromArgb(0, 200, 0));
            g.AddGrip(0.75, Color.FromArgb(220, 220, 0));
            g.AddGrip(1.00, Color.FromArgb(220, 0, 0));
            return g;
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("import.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("3C7F1E08-5B62-4A9D-91C4-7E2A60D5F3B1"); }
        }
    }
}
