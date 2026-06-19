using System;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Rhino.Geometry;

namespace bertini_real
{
    /// <summary>
    /// Reads a self-contained standalone-curve export (br_gh_export.json) written by Python's
    /// Curve.export_gh_json.  Brings the vertices in as ONE unified set and exposes each curve
    /// piece as a polyline referring to those vertices by index.
    /// </summary>
    public class CurveReadGhJson : GH_Component
    {
        public CurveReadGhJson()
          : base("Curve Read GH JSON", "CurveReadJSON",
                 "Read a bertini_real curve export: one unified vertex set, curve pieces as polylines",
                 "bertini_real", "Curve")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("File Path", "F", "Path to br_gh_export.json (a curve export)", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Vertices", "V", "The single unified vertex set; list index = global vertex id", GH_ParamAccess.list);
            pManager.AddCurveParameter("Curves", "C", "One polyline per curve piece", GH_ParamAccess.tree);
            pManager.AddTextParameter("Curve Types", "T", "Type tag per curve piece, parallel to Curves", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("Curve Indices", "CI", "Per curve piece: vertex indices into Vertices", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            string path = "";
            if (!DA.GetData(0, ref path)) return;

            GhExport content;
            try
            {
                content = GhJsonIO.Load(path);
            }
            catch (Exception e)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Failed to read/parse JSON: " + e.Message);
                return;
            }

            if (content == null || content.decomposition_type != "curve")
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Not a curve export (decomposition_type != 'curve').");
                return;
            }

            var verts = GhJsonIO.ToVertices(content);

            var curves = new DataTree<Curve>();
            var types = new DataTree<string>();
            var curveIdx = new DataTree<int>();

            if (content.curve_pieces != null)
            {
                foreach (var cp in content.curve_pieces)
                {
                    var branch = new GH_Path(cp.piece_index);

                    PolylineCurve pl = GhJsonIO.ToPolyline(cp.vertex_indices, verts);
                    if (pl != null) curves.Add(pl, branch);

                    types.Add(cp.type, branch);
                    if (cp.vertex_indices != null) curveIdx.AddRange(cp.vertex_indices, branch);
                }
            }

            DA.SetDataList(0, verts);
            DA.SetDataTree(1, curves);
            DA.SetDataTree(2, types);
            DA.SetDataTree(3, curveIdx);
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("telephone.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("7C2A9F10-3B8E-4D55-A1C7-6E0B4F92D38B"); }
        }
    }
}
