using System;
using System.Collections.Generic;
using System.IO;
using System.Text.Json;
using Rhino.Geometry;

namespace bertini_real
{
    /// <summary>
    /// Shared, Rhino-light helpers for reading the self-contained br_gh_export.json and turning
    /// its unified vertex set + index data into Rhino geometry.  Kept separate from the GH
    /// components so the mapping logic stays small and could be unit-tested later.
    /// </summary>
    internal static class GhJsonIO
    {
        public static GhExport Load(string path)
        {
            return JsonSerializer.Deserialize<GhExport>(File.ReadAllText(path));
        }

        /// <summary>The single unified vertex set, as Rhino points (index = global vertex id).</summary>
        public static List<Point3d> ToVertices(GhExport content)
        {
            var verts = new List<Point3d>();
            if (content?.vertices == null) return verts;
            foreach (var v in content.vertices)
            {
                double x = v != null && v.Length > 0 ? v[0] : 0.0;
                double y = v != null && v.Length > 1 ? v[1] : 0.0;
                double z = v != null && v.Length > 2 ? v[2] : 0.0;
                verts.Add(new Point3d(x, y, z));
            }
            return verts;
        }

        /// <summary>
        /// Build a mesh on the FULL shared vertex cloud, so triangle indices stay global and
        /// coincident vertices across pieces are recognized as the same point (enabling joins,
        /// closed-solid detection, exact curve/mesh intersections).  Deliberately does NOT
        /// cull/compact vertices -- that would reindex and destroy the shared identity.
        /// Degenerate and out-of-range triangles are skipped without dropping vertices.
        /// </summary>
        public static Mesh BuildMesh(GhMesh g, List<Point3d> verts)
        {
            if (g?.triangles == null) return null;

            var mesh = new Mesh();
            mesh.Vertices.AddVertices(verts);

            int n = verts.Count;
            int[] tri = g.triangles;
            for (int t = 0; t + 2 < tri.Length; t += 3)
            {
                int a = tri[t], b = tri[t + 1], c = tri[t + 2];
                if (a < 0 || b < 0 || c < 0 || a >= n || b >= n || c >= n) continue; // out of range
                if (a == b || b == c || a == c) continue;                            // degenerate
                mesh.Faces.AddFace(a, b, c);
            }

            mesh.Normals.ComputeNormals();
            return mesh;
        }

        /// <summary>
        /// Polyline through the shared vertices, by index.  Returns null for &lt;2 valid points
        /// (e.g. a nodal singularity), which the caller skips.
        /// </summary>
        public static PolylineCurve ToPolyline(int[] indices, List<Point3d> verts)
        {
            if (indices == null) return null;

            var pts = new List<Point3d>();
            int n = verts.Count;
            foreach (int i in indices)
                if (i >= 0 && i < n) pts.Add(verts[i]);

            if (pts.Count < 2) return null;
            return new PolylineCurve(pts);
        }

        /// <summary>
        /// The decomposition's bounding sphere as a closed Brep (ready for boolean / capping
        /// operations), or null if the sphere data is missing or degenerate.
        /// </summary>
        public static Brep ToSphereBrep(GhSphere s)
        {
            if (s == null || s.center == null || s.radius <= 0.0) return null;

            double[] c = s.center;
            var center = new Point3d(
                c.Length > 0 ? c[0] : 0.0,
                c.Length > 1 ? c[1] : 0.0,
                c.Length > 2 ? c[2] : 0.0);

            return new Sphere(center, s.radius).ToBrep();
        }
    }
}
