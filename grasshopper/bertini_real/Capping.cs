using System;
using System.Collections.Generic;
using Rhino;
using Rhino.Geometry;

namespace bertini_real
{
    /// <summary>
    /// Shared helpers for the cap components (Sphere Caps, Flat Caps): finding the piece mesh's
    /// own naked boundary loops that lie on the sphere, and adding non-degenerate triangles.
    /// </summary>
    internal static class Capping
    {
        // a tiny absolute distance below which two cap vertices are treated as the same point
        public const double CoincidentTol = 1e-9;

        /// <summary>
        /// Closed loops (ordered lists of TopologyVertex indices) of the mesh's naked boundary that
        /// lie on the sphere.  Sets <paramref name="unclean"/> if some on-sphere boundary did not
        /// form clean degree-2 cycles.  Works in topology space (coincident vertices merged); zero-
        /// length edges collapse and are ignored.
        /// </summary>
        public static List<List<int>> OnSphereBoundaryLoops(Mesh mesh, Point3d center, double radius, double tol, out bool unclean)
        {
            unclean = false;
            var result = new List<List<int>>();

            var topo = mesh.TopologyVertices;
            var edges = mesh.TopologyEdges;

            bool OnSphere(int tv)
            {
                Point3d p = topo[tv];
                return Math.Abs(p.DistanceTo(center) - radius) < tol;
            }

            var adj = new Dictionary<int, List<int>>();
            void Link(int u, int v)
            {
                if (!adj.TryGetValue(u, out var lu)) { lu = new List<int>(); adj[u] = lu; }
                if (!lu.Contains(v)) lu.Add(v);
            }

            for (int e = 0; e < edges.Count; e++)
            {
                if (edges.GetConnectedFaces(e).Length != 1) continue; // not naked
                IndexPair ip = edges.GetTopologyVertices(e);
                int i = ip.I, j = ip.J;
                if (i == j) continue;                 // degenerate (collapsed) edge
                if (!OnSphere(i) || !OnSphere(j)) continue;
                Link(i, j);
                Link(j, i);
            }

            if (adj.Count == 0) return result;

            var visited = new HashSet<int>();
            foreach (int start in adj.Keys)
            {
                if (visited.Contains(start)) continue;

                var loop = new List<int>();
                int prev = -1, cur = start;
                bool clean = true;

                while (true)
                {
                    visited.Add(cur);
                    loop.Add(cur);

                    var nbrs = adj[cur];
                    if (nbrs.Count != 2) { clean = false; break; } // junction or dead-end

                    int next = nbrs[0] != prev ? nbrs[0] : nbrs[1];
                    if (next == start) break;          // closed the loop
                    if (visited.Contains(next)) { clean = false; break; }
                    prev = cur;
                    cur = next;
                }

                if (clean && loop.Count >= 3)
                    result.Add(loop);
                else
                    unclean = true;
            }

            return result;
        }

        /// <summary>
        /// Add a triangle, skipping only TRULY degenerate ones (a collapsed edge -- two coincident
        /// vertices); thin-but-valid triangles must be kept, or high-resolution caps develop holes.
        /// </summary>
        public static void AddTri(Mesh m, int i, int j, int k)
        {
            Point3d a = m.Vertices[i];
            Point3d b = m.Vertices[j];
            Point3d c = m.Vertices[k];
            if (a.DistanceTo(b) < CoincidentTol ||
                b.DistanceTo(c) < CoincidentTol ||
                a.DistanceTo(c) < CoincidentTol)
                return;
            m.Faces.AddFace(i, j, k);
        }
    }
}
