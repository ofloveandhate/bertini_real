using System;
using System.Collections.Generic;
using System.Linq;
using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

namespace bertini_real
{
    /// <summary>
    /// Spreads pieces along their connector axes to show how they would assemble -- a directional
    /// assembly explosion, as opposed to the radial Spread Pieces.
    ///
    /// The pieces and their singularities form a graph (each singularity is an edge joining exactly
    /// two pieces, with a connector direction).  This does a breadth-first traversal from a fixed
    /// root: each piece's displacement is its parent's displacement PLUS a step of Distance along
    /// the connecting axis, oriented so the child slides off the rod away from the parent.  So a
    /// chain telescopes outward and the root stays put -- robust for trees and chains, and it
    /// spanning-trees any cycles.  The root defaults to the most-connected piece (a natural hub);
    /// set Root to fix a specific piece.  Disconnected groups are each rooted independently.
    ///
    /// Wire Sing Locations / Sing Directions / Sing On Pieces straight from the reader.
    /// </summary>
    public class SurfaceSpreadByConnectors : GH_Component
    {
        public SurfaceSpreadByConnectors()
          : base("Spread By Connectors", "SpreadConn",
                 "Assembly explosion: BFS the connector graph and slide each piece off its rods along the connector axes",
                 "bertini_real", "Surface")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGeometryParameter("Geometry", "G", "Per-piece geometry, one branch per piece (meshes, or mesh + connectors)", GH_ParamAccess.tree);
            pManager.AddPointParameter("Sing Locations", "SL", "Singularity locations (from Surface Read GH JSON)", GH_ParamAccess.list);
            pManager.AddVectorParameter("Sing Directions", "SD", "Singularity connector directions (from Surface Read GH JSON)", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Sing On Pieces", "SOP", "Per piece: indices of the singularities on it (from Surface Read GH JSON)", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Distance", "D", "Separation distance per connector along the assembly chain (model units). 0 = no move.", GH_ParamAccess.item, 1.0);
            pManager.AddIntegerParameter("Root", "R", "Piece index to hold fixed (-1 = auto: the most-connected piece)", GH_ParamAccess.item, -1);
            Params.Input[4].Optional = true;
            Params.Input[5].Optional = true;
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGeometryParameter("Geometry", "G", "Spread-apart geometry (same tree structure)", GH_ParamAccess.tree);
            pManager.AddVectorParameter("Translations", "T", "Translation applied to each piece", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<IGH_GeometricGoo> geometry;
            if (!DA.GetDataTree(0, out geometry)) return;

            var locations = new List<Point3d>();
            var directions = new List<Vector3d>();
            GH_Structure<GH_Integer> onPieces;
            DA.GetDataList(1, locations);
            DA.GetDataList(2, directions);
            if (!DA.GetDataTree(3, out onPieces)) return;

            double distance = 1.0;
            DA.GetData(4, ref distance);
            int root = -1;
            DA.GetData(5, ref root);

            // pieces we have geometry for, with their centers and branch paths
            var pieces = new List<int>();
            var centroid = new Dictionary<int, Point3d>();
            var branchPath = new Dictionary<int, GH_Path>();
            for (int b = 0; b < geometry.PathCount; b++)
            {
                GH_Path path = geometry.get_Path(b);
                int pi = path.Indices.Length > 0 ? path.Indices[path.Indices.Length - 1] : b;
                BoundingBox bb = BoundingBox.Empty;
                bool any = false;
                foreach (var goo in geometry.get_Branch(path))
                    if (goo is IGH_GeometricGoo gg && gg.IsValid) { bb.Union(gg.Boundingbox); any = true; }
                if (!any || centroid.ContainsKey(pi)) continue;
                pieces.Add(pi);
                centroid[pi] = bb.Center;
                branchPath[pi] = path;
            }
            if (pieces.Count == 0) return;

            // invert Sing On Pieces -> which pieces each singularity touches (restricted to our pieces)
            var singToPieces = new Dictionary<int, List<int>>();
            foreach (int pi in pieces)
            {
                var sopPath = new GH_Path(pi);
                if (!onPieces.PathExists(sopPath)) continue;
                foreach (var goo in onPieces.get_Branch(sopPath))
                    if (goo is GH_Integer gi)
                    {
                        if (!singToPieces.TryGetValue(gi.Value, out var lst)) { lst = new List<int>(); singToPieces[gi.Value] = lst; }
                        if (!lst.Contains(pi)) lst.Add(pi);
                    }
            }

            // adjacency: a singularity joining exactly two pieces is an edge
            var adj = new Dictionary<int, List<(int nbr, int sing)>>();
            void Link(int a, int bb2, int s)
            {
                if (!adj.TryGetValue(a, out var la)) { la = new List<(int, int)>(); adj[a] = la; }
                la.Add((bb2, s));
            }
            foreach (var kv in singToPieces)
                if (kv.Value.Count == 2)
                {
                    Link(kv.Value[0], kv.Value[1], kv.Key);
                    Link(kv.Value[1], kv.Value[0], kv.Key);
                }

            // BFS from a root, accumulating displacement along the connector chain
            var disp = new Dictionary<int, Vector3d>();
            foreach (int pi in pieces) disp[pi] = Vector3d.Zero;
            var visited = new HashSet<int>();

            // root order: an explicit Root first, then most-connected pieces (one root per component)
            var order = pieces.OrderByDescending(pi => adj.TryGetValue(pi, out var l) ? l.Count : 0).ToList();
            if (root >= 0 && centroid.ContainsKey(root)) { order.Remove(root); order.Insert(0, root); }

            int components = 0;
            foreach (int start in order)
            {
                if (visited.Contains(start)) continue;
                components++;
                var queue = new Queue<int>();
                queue.Enqueue(start);
                visited.Add(start);
                disp[start] = Vector3d.Zero;

                while (queue.Count > 0)
                {
                    int cur = queue.Dequeue();
                    if (!adj.TryGetValue(cur, out var nbrs)) continue;
                    foreach (var (nbr, sing) in nbrs)
                    {
                        if (visited.Contains(nbr)) continue;
                        if (sing < 0 || sing >= directions.Count || sing >= locations.Count) continue;

                        Vector3d axis = directions[sing];
                        if (axis.IsTiny()) continue;
                        axis.Unitize();

                        // orient the axis so the child slides away from the singularity (and parent)
                        double along = (centroid[nbr] - locations[sing]) * axis;
                        Vector3d step = distance * (along >= 0 ? 1.0 : -1.0) * axis;

                        disp[nbr] = disp[cur] + step;
                        visited.Add(nbr);
                        queue.Enqueue(nbr);
                    }
                }
            }

            if (components > 1)
                AddRuntimeMessage(GH_RuntimeMessageLevel.Remark,
                    $"{components} disconnected connector groups; each is rooted independently and may overlap at the origin.");

            // emit
            var outGeo = new DataTree<IGH_GeometricGoo>();
            var outVec = new DataTree<Vector3d>();
            foreach (int pi in pieces)
            {
                GH_Path path = branchPath[pi];
                Transform xf = Transform.Translation(disp[pi]);
                foreach (var goo in geometry.get_Branch(path))
                {
                    if (!(goo is IGH_GeometricGoo gg) || !gg.IsValid) continue;
                    IGH_GeometricGoo moved = gg.DuplicateGeometry();
                    moved = moved.Transform(xf);
                    outGeo.Add(moved, path);
                }
                outVec.Add(disp[pi], path);
            }

            DA.SetDataTree(0, outGeo);
            DA.SetDataTree(1, outVec);
        }

        protected override System.Drawing.Bitmap Icon => IconLoader.GetIcon("transform.png");

        public override Guid ComponentGuid
        {
            get { return new Guid("9F26C7B4-3A81-4D60-8E15-7C0B2F95E6A3"); }
        }
    }
}
