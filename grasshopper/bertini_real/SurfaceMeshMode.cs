using System;
using Grasshopper.Kernel;

namespace bertini_real
{
    /// <summary>
    /// Maps an integer to a Surface Read GH JSON "Mesh Mode" string, so a numeric slider can pick
    /// the mode: 0 = auto, 1 = smooth, 2 = raw.  Feed the output into the reader's Mesh Mode input.
    /// </summary>
    public class SurfaceMeshMode : GH_Component
    {
        private static readonly string[] Modes = { "auto", "smooth", "raw" };

        public SurfaceMeshMode()
          : base("Mesh Mode", "MeshMode",
                 "Pick a Surface Read GH JSON Mesh Mode by index: 0 = auto, 1 = smooth, 2 = raw",
                 "bertini_real", "Surface")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddIntegerParameter("Index", "i", "0 = auto, 1 = smooth, 2 = raw", GH_ParamAccess.item, 0);
            Params.Input[0].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Mesh Mode", "MM", "Mode string for Surface Read GH JSON", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            int index = 0;
            DA.GetData(0, ref index);

            int clamped = Math.Max(0, Math.Min(Modes.Length - 1, index));
            if (clamped != index)
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                    $"Index {index} out of range [0..{Modes.Length - 1}]; using '{Modes[clamped]}'.");

            DA.SetData(0, Modes[clamped]);
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("import.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("2E5D8B30-7A19-4C62-9F41-0B6C3E8A75D2"); }
        }
    }
}
