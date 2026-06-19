using System;
using System.Collections.Generic;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Rhino.Geometry;

namespace bertini_real
{
    /// <summary>
    /// Reads a self-contained surface export (br_gh_export.json) written by Python's
    /// Surface.export_gh_json.  Brings the vertices in as ONE unified set, then exposes each
    /// nonsingular piece as a mesh and the curve pieces embedded on it as polylines -- all
    /// referring to the same vertices by index.
    /// </summary>
    public class SurfaceReadGhJson : GH_Component
    {
        public SurfaceReadGhJson()
          : base("Surface Read GH JSON", "SurfReadJSON",
                 "Read a bertini_real surface export: one unified vertex set, pieces as meshes, embedded curves as polylines",
                 "bertini_real", "Surface")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("File Path", "F", "Path to br_gh_export.json (a surface export)", GH_ParamAccess.item);
            pManager.AddTextParameter("Mesh Mode", "MM", "auto | smooth | raw  (auto = smooth when sampled, else raw)", GH_ParamAccess.item, "auto");
            Params.Input[1].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Vertices", "V", "The single unified vertex set; list index = global vertex id", GH_ParamAccess.list);
            pManager.AddMeshParameter("Meshes", "M", "One mesh per nonsingular piece, built on the shared vertices", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("Mesh Faces", "MF", "Per piece: flat triangle indices into Vertices", GH_ParamAccess.tree);
            pManager.AddCurveParameter("Curves", "C", "Embedded curve pieces per surface piece", GH_ParamAccess.tree);
            pManager.AddTextParameter("Curve Types", "T", "Type tag per curve (critical/sphere/singular/midslice/critslice), parallel to Curves", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("Curve Indices", "CI", "Per curve: vertex indices into Vertices (path {piece, curve})", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("Face Indices", "FI", "Global surface face ids per piece", GH_ParamAccess.tree);
            pManager.AddBrepParameter("Sphere", "S", "Bounding sphere of the decomposition as a closed Brep", GH_ParamAccess.item);
            pManager.AddPointParameter("Sing Locations", "SL", "Nodal singularity locations (one per singularity)", GH_ParamAccess.list);
            pManager.AddVectorParameter("Sing Directions", "SD", "Nodal singularity connector axis directions (one per singularity)", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Sing Parities", "SP", "Per singularity: parity (-1/0/1) on each piece", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("Sing On Pieces", "SOP", "Per piece: indices of the singularities on it", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            string path = "";
            string mode = "auto";
            if (!DA.GetData(0, ref path)) return;
            DA.GetData(1, ref mode);
            mode = (mode ?? "auto").Trim().ToLowerInvariant();

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

            if (content == null || content.decomposition_type != "surface")
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Not a surface export (decomposition_type != 'surface').");
                return;
            }

            var verts = GhJsonIO.ToVertices(content);

            var meshes = new DataTree<Mesh>();
            var meshFaces = new DataTree<int>();
            var curves = new DataTree<Curve>();
            var types = new DataTree<string>();
            var curveIdx = new DataTree<int>();
            var faceIdx = new DataTree<int>();

            if (content.pieces != null)
            {
                foreach (var piece in content.pieces)
                {
                    var branch = new GH_Path(piece.piece_index);

                    GhMesh chosen = PickMesh(piece, mode);
                    if (chosen?.triangles != null)
                    {
                        Mesh m = GhJsonIO.BuildMesh(chosen, verts);
                        if (m != null) meshes.Add(m, branch);
                        meshFaces.AddRange(chosen.triangles, branch);
                    }
                    else
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, $"Piece {piece.piece_index} has no usable mesh.");
                    }

                    if (piece.face_indices != null)
                        faceIdx.AddRange(piece.face_indices, branch);

                    if (piece.curves != null)
                    {
                        int ord = 0;
                        foreach (var c in piece.curves)
                        {
                            PolylineCurve pl = GhJsonIO.ToPolyline(c.vertex_indices, verts);
                            if (pl == null) continue; // <2 points (e.g. nodal singularity)

                            curves.Add(pl, branch);
                            types.Add(c.type, branch);
                            curveIdx.AddRange(c.vertex_indices ?? Array.Empty<int>(), new GH_Path(piece.piece_index, ord));
                            ord++;
                        }
                    }
                }
            }

            DA.SetDataList(0, verts);
            DA.SetDataTree(1, meshes);
            DA.SetDataTree(2, meshFaces);
            DA.SetDataTree(3, curves);
            DA.SetDataTree(4, types);
            DA.SetDataTree(5, curveIdx);
            DA.SetDataTree(6, faceIdx);

            Brep sphere = GhJsonIO.ToSphereBrep(content.sphere);
            if (sphere != null)
                DA.SetData(7, sphere);
            else
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No valid bounding sphere in the export.");

            // singularity / connector data (folded into the same JSON)
            var singLocations = new List<Point3d>();
            var singDirections = new List<Vector3d>();
            var singParities = new DataTree<int>();
            var singOnPieces = new DataTree<int>();

            var sg = content.singularities;
            if (sg != null)
            {
                if (sg.locations != null)
                    foreach (var p in sg.locations)
                        if (p != null && p.Length >= 3)
                            singLocations.Add(new Point3d(p[0], p[1], p[2]));

                if (sg.directions != null)
                    foreach (var d in sg.directions)
                        if (d != null && d.Length >= 3)
                            singDirections.Add(new Vector3d(d[0], d[1], d[2]));

                if (sg.parities != null)
                    for (int s = 0; s < sg.parities.Length; s++)
                        if (sg.parities[s] != null)
                            singParities.AddRange(sg.parities[s], new GH_Path(s));

                if (sg.on_pieces != null)
                    for (int pc = 0; pc < sg.on_pieces.Length; pc++)
                        if (sg.on_pieces[pc] != null)
                            singOnPieces.AddRange(sg.on_pieces[pc], new GH_Path(pc));
            }

            DA.SetDataList(8, singLocations);
            DA.SetDataList(9, singDirections);
            DA.SetDataTree(10, singParities);
            DA.SetDataTree(11, singOnPieces);
        }

        private static GhMesh PickMesh(GhPiece piece, string mode)
        {
            if (mode == "raw") return piece.mesh_raw;
            if (mode == "smooth") return piece.mesh_smooth ?? piece.mesh_raw;
            return piece.mesh_smooth ?? piece.mesh_raw; // auto
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("lego.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("B5E3C7A2-1D4F-4E8A-9C6B-2F7A0D9E13A4"); }
        }
    }
}
