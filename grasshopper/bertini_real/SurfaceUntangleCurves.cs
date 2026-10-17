using System;
using System.Collections;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

namespace bertini_real
{
    /// <summary>
    /// Splits the parallel Curves + Curve Types trees coming out of "Surface Read GH JSON" into
    /// one output per curve type, so you can grab just the type you want without manual tree
    /// filtering / get-item juggling.  Per-piece tree paths are preserved on every output, so a
    /// curve stays associated with the piece it came from.
    /// </summary>
    public class SurfaceUntangleCurves : GH_Component
    {
        public SurfaceUntangleCurves()
          : base("Untangle Curves By Type", "UntangleCurves",
                 "Split the embedded curves from Surface Read GH JSON into one output per type",
                 "bertini_real", "Surface")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Curves", "C", "Embedded curves tree from Surface Read GH JSON", GH_ParamAccess.tree);
            pManager.AddTextParameter("Curve Types", "T", "Curve Types tree from Surface Read GH JSON (parallel to Curves)", GH_ParamAccess.tree);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddCurveParameter("Critical", "Cr", "Critical curve pieces", GH_ParamAccess.tree);
            pManager.AddCurveParameter("Sphere", "Sp", "Sphere curve pieces", GH_ParamAccess.tree);
            pManager.AddCurveParameter("Singular", "Si", "Singular curve pieces", GH_ParamAccess.tree);
            pManager.AddCurveParameter("Midslice", "Mi", "Midslice curve pieces", GH_ParamAccess.tree);
            pManager.AddCurveParameter("Critslice", "Cs", "Critslice curve pieces", GH_ParamAccess.tree);
            pManager.AddCurveParameter("Other", "Ot", "Curves whose type is unknown / unrecognized", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<GH_Curve> curves;
            GH_Structure<GH_String> types;
            if (!DA.GetDataTree(0, out curves)) return;
            if (!DA.GetDataTree(1, out types)) return;

            var critical = new GH_Structure<GH_Curve>();
            var sphere = new GH_Structure<GH_Curve>();
            var singular = new GH_Structure<GH_Curve>();
            var midslice = new GH_Structure<GH_Curve>();
            var critslice = new GH_Structure<GH_Curve>();
            var other = new GH_Structure<GH_Curve>();

            for (int b = 0; b < curves.PathCount; b++)
            {
                GH_Path path = curves.get_Path(b);
                IList curveBranch = curves.get_Branch(path);
                IList typeBranch = types.PathExists(path) ? types.get_Branch(path) : null;

                if (typeBranch == null || typeBranch.Count != curveBranch.Count)
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                        $"Types do not line up with Curves at path {path}; missing tags treated as 'Other'.");

                for (int i = 0; i < curveBranch.Count; i++)
                {
                    var gc = curveBranch[i] as GH_Curve;
                    if (gc == null) continue;

                    string t = "unknown";
                    if (typeBranch != null && i < typeBranch.Count && typeBranch[i] is GH_String gs && gs.Value != null)
                        t = gs.Value.Trim().ToLowerInvariant();

                    GH_Structure<GH_Curve> target;
                    switch (t)
                    {
                        case "critical": target = critical; break;
                        case "sphere": target = sphere; break;
                        case "singular": target = singular; break;
                        case "midslice": target = midslice; break;
                        case "critslice": target = critslice; break;
                        default: target = other; break;
                    }

                    target.Append(gc, path);
                }
            }

            DA.SetDataTree(0, critical);
            DA.SetDataTree(1, sphere);
            DA.SetDataTree(2, singular);
            DA.SetDataTree(3, midslice);
            DA.SetDataTree(4, critslice);
            DA.SetDataTree(5, other);
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("transform.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("3F1B8D6C-9A42-4C71-B8E5-71D0A6F4C982"); }
        }
    }
}
