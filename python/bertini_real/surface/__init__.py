"""
    :platform: Unix, Windows
    :synopsis: This module contains Surface and SurfacePiece types.

"""

# Foong Min Wong
# Fall 2018 - Spring 2019
#
# Silviana Amethyst
# Spring 2022, Summer 2022
#
# University of Wisconsin, Eau Claire

import bertini_real.parse
import bertini_real.exception as br_except
import numpy as np
from bertini_real.decomposition import Decomposition
from bertini_real.curve import Curve, CurvePiece, is_edge_degenerate, _points_to_xyz
from bertini_real.vertex import Vertex
from bertini_real.vertex import VertexType
from bertini_real.util import ReversableList

import os
import enum
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.widgets import CheckButtons

from collections import defaultdict

import copy
import math
import trimesh

import json





_default_file_type = 'stl'

_default_solidify_thickness =  0.1

_default_piece_basename_smooth = 'br_piece_smooth'
_default_piece_basename_raw = 'br_piece_raw'

_default_surface_basename_smooth = 'br_surface_smooth'
_default_surface_basename_raw = 'br_surface_raw'






def _mesh_triangles(mesh):
    """
    flatten a `trimesh.Trimesh`'s faces into a dict for the Grasshopper JSON export.

    the triangle entries are indices into the surface's unified (global) vertex set,
    because `as_mesh_raw`/`as_mesh_smooth` build the mesh from `extract_points()` with
    `keep_all_vertices=True` (so trimesh does not reindex).  returns None if mesh is None.
    """
    if mesh is None:
        return None

    triangles = np.asarray(mesh.faces, dtype=int).reshape(-1).tolist()
    return {
        "triangles": triangles,
        "triangle_count": len(mesh.faces),
    }


def export_mesh(mesh, basename, autoname_using_folder=False, file_type=_default_file_type, verbose=True):
    """
    Saves a mesh (generated elsewhere) to disk,
    and returns the name of the file which was saved
    """


    if autoname_using_folder:
        fileName = os.getcwd().split(os.sep)[-1]
        outname = f'{basename}_{fileName}.{file_type}'
    else:
        outname = f'{basename}.{file_type}'

    mesh.export(file_obj=outname, file_type=file_type)

    if verbose:
        print("Exported \x1b[0;35;40m " + outname + "\x1b[0m successfully")

    return outname


def solidify_mesh(mesh, distance, offset=0):

    A = mesh # relying on referenced nature of Python here

    A.fix_normals()

    B = copy.deepcopy(mesh)

    # reverse every triangles and flip every normals
    B.invert()

    # calculate A, B vertex normals
    vertexnormsA = A.vertex_normals
    vertexnormsB = B.vertex_normals

    distA = (distance) * (offset + 1) / 2
    distB = (distance) * (1 - (offset + 1) / 2)

    # create A & B vertices that move corresponding to vertex normals and
    # distance
    A.vertices = [v + vn * distA for v,
                  vn in zip(A.vertices, A.vertex_normals)]
    B.vertices = [v + vn * distB for v,
                  vn in zip(B.vertices, B.vertex_normals)]

    numVerts = len(A.vertices)

    boundary_triangles = []

    boundary_groups = trimesh.grouping.group_rows(
        A.edges_sorted, require_count=1)

    boundary_edges = A.edges[boundary_groups]

    for edge in boundary_edges:
        for i in range(len(edge) - 1):
            t1 = [edge[i+1], edge[i], edge[i] + numVerts]
            t2 = [edge[i + 1], edge[i] + numVerts, edge[i + 1] + numVerts]

            boundary_triangles.append(t1)
            boundary_triangles.append(t2)

    Q = np.concatenate((A.vertices, B.vertices), axis=0)

    newBoundary = trimesh.Trimesh(Q, boundary_triangles)

    finalmesh = A + newBoundary + B
    return finalmesh


# i used the following to help me solve this problem:
#     https://stackoverflow.com/questions/61531935/irerate-over-package-data-files-and-copy-them-to-current-working-directory
#    lol, the misspelled "iterate" is not my fault.
#
# also:
#     https://stackoverflow.com/questions/32490629/getting-todays-date-in-yyyy-mm-dd-in-python
def copy_all_scad_files_here():
    """
    copy all source .scad files provided in bertini_real to the current directory
    """

    import pkgutil
    import pkg_resources
    from os.path import join

    scad_files = pkg_resources.resource_listdir("bertini_real", "surface/scad")

    for s in scad_files:
        contents = pkgutil.get_data('bertini_real',join('surface/scad/',s))

        with open(s,'wb') as f:
            f.write(contents)




def _slerp(v0, v1, t):
    """Spherical interpolation of two unit vectors; linear fallback when (anti)parallel."""
    dot = np.clip(np.dot(v0, v1), -1.0, 1.0)
    omega = np.arccos(dot)
    so = np.sin(omega)
    if omega < 1e-9 or so < 1e-9:
        lin = (1.0 - t) * v0 + t * v1
        n = np.linalg.norm(lin)
        return v0 if n < 1e-12 else lin / n
    return np.sin((1.0 - t) * omega) / so * v0 + np.sin(t * omega) / so * v1


def _on_sphere_boundary_loops(mesh, center, radius, tol):
    """
    Ordered loops (lists of vertex indices) of `mesh`'s naked boundary edges that lie on the
    sphere of the given center/radius.  Only clean degree-2 cycles are returned; open arcs or
    non-manifold junctions are dropped.
    """
    from collections import Counter, defaultdict

    V = np.asarray(mesh.vertices)
    edge_count = Counter()
    for tri in mesh.faces:
        a, b, c = int(tri[0]), int(tri[1]), int(tri[2])
        for u, v in ((a, b), (b, c), (c, a)):
            edge_count[frozenset((u, v))] += 1
    naked = [tuple(e) for e, cnt in edge_count.items() if cnt == 1 and len(e) == 2]

    def on_sphere(i):
        return abs(np.linalg.norm(V[i][:3] - center) - radius) < tol

    adj = defaultdict(list)
    for e in naked:
        a, b = tuple(e)
        if a != b and on_sphere(a) and on_sphere(b):
            adj[a].append(b)
            adj[b].append(a)

    loops = []
    seen = set()
    for start in list(adj):
        if start in seen:
            continue
        loop = [start]
        seen.add(start)
        prev, cur, ok = -1, start, True
        while True:
            nbrs = adj[cur]
            if len(nbrs) != 2:
                ok = False
                break
            nxt = nbrs[0] if nbrs[0] != prev else nbrs[1]
            if nxt == start:
                break
            if nxt in seen:
                ok = False
                break
            loop.append(nxt)
            seen.add(nxt)
            prev, cur = cur, nxt
        if ok and len(loop) >= 3:
            loops.append(loop)
    return loops


def sphere_cap_meshes(mesh, center, radius, resolution=4, tol=1e-3):
    """
    Faceted spherical caps that close `mesh`'s naked boundary loops lying on the sphere.

    The Python twin of the Grasshopper "Sphere Caps" component: it caps the mesh's OWN
    boundary (so the caps share its vertices and weld watertight), keeps the smaller-area
    side of each loop, and subdivides each cap into `resolution` radial rings slerped along
    the sphere.  Returns a list of `trimesh.Trimesh`.
    """
    center = np.asarray(center, dtype=float)[:3]
    V = np.asarray(mesh.vertices)
    caps = []

    for loop in _on_sphere_boundary_loops(mesh, center, radius, tol):
        pts = np.array([V[i][:3] for i in loop])
        dirs = np.array([(p - center) / np.linalg.norm(p - center) for p in pts])
        mean = pts.mean(axis=0) - center
        nrm = np.linalg.norm(mean)
        mean = mean / nrm if nrm > 1e-12 else np.array([0.0, 0.0, 1.0])

        def fan_area(sign):
            apex = center + sign * radius * mean
            return sum(np.linalg.norm(np.cross(pts[k] - apex, pts[(k + 1) % len(pts)] - apex)) / 2.0
                       for k in range(len(pts)))

        sign = 1 if fan_area(1) <= fan_area(-1) else -1
        pole_dir = sign * mean
        pole = center + radius * pole_dir

        R = max(1, int(resolution))
        n = len(loop)
        verts = [p for p in pts]
        for r in range(1, R):
            t = r / R
            for k in range(n):
                verts.append(center + radius * _slerp(dirs[k], pole_dir, t))
        pole_i = len(verts)
        verts.append(pole)
        verts = np.array(verts)

        def vid(r, k):
            return r * n + k

        faces = []

        def add(i, j, k):
            # skip only truly degenerate (collapsed-edge) triangles; keep thin ones
            if (np.linalg.norm(verts[i] - verts[j]) < 1e-9 or
                    np.linalg.norm(verts[j] - verts[k]) < 1e-9 or
                    np.linalg.norm(verts[i] - verts[k]) < 1e-9):
                return
            faces.append([i, j, k])

        for r in range(R - 1):
            for k in range(n):
                k2 = (k + 1) % n
                add(vid(r, k), vid(r, k2), vid(r + 1, k2))
                add(vid(r, k), vid(r + 1, k2), vid(r + 1, k))
        for k in range(n):
            add(vid(R - 1, k), vid(R - 1, (k + 1) % n), pole_i)

        if faces:
            caps.append(trimesh.Trimesh(verts, np.array(faces), process=False))

    return caps


def flat_cap_meshes(mesh, center, radius, resolution=1, tol=1e-3):
    """
    Flat caps that close `mesh`'s on-sphere boundary loops with a fan to each loop's centroid --
    the flat alternative to sphere_cap_meshes, and the Python twin of the "Flat Caps" component.

    Same on-sphere boundary as sphere_cap_meshes, but the apex is the loop centroid (a flat fill)
    rather than a point on the sphere, so e.g. a Bertini cylinder gets flat disk ends.  Caps the
    mesh's OWN boundary (welds watertight), with `resolution` concentric linearly-interpolated
    rings (1 = a single flat fan).  Returns a list of trimesh.Trimesh.
    """
    center = np.asarray(center, dtype=float)[:3]
    V = np.asarray(mesh.vertices)
    caps = []

    for loop in _on_sphere_boundary_loops(mesh, center, radius, tol):
        pts = np.array([V[i][:3] for i in loop])
        centroid = pts.mean(axis=0)

        R = max(1, int(resolution))
        n = len(loop)
        verts = [p for p in pts]
        for r in range(1, R):
            t = r / R
            for k in range(n):
                verts.append(pts[k] + t * (centroid - pts[k]))
        apex_i = len(verts)
        verts.append(centroid)
        verts = np.array(verts)

        def vid(r, k):
            return r * n + k

        faces = []

        def add(i, j, k):
            if (np.linalg.norm(verts[i] - verts[j]) < 1e-9 or
                    np.linalg.norm(verts[j] - verts[k]) < 1e-9 or
                    np.linalg.norm(verts[i] - verts[k]) < 1e-9):
                return
            faces.append([i, j, k])

        for r in range(R - 1):
            for k in range(n):
                k2 = (k + 1) % n
                add(vid(r, k), vid(r, k2), vid(r + 1, k2))
                add(vid(r, k), vid(r + 1, k2), vid(r + 1, k))
        for k in range(n):
            add(vid(R - 1, k), vid(R - 1, (k + 1) % n), apex_i)

        if faces:
            caps.append(trimesh.Trimesh(verts, np.array(faces), process=False))

    return caps


def join_meshes(meshes):
    """
    Concatenate meshes and merge coincident vertices into one (ideally watertight) trimesh.
    The Python twin of the Grasshopper "Close Piece" weld.  Returns None if nothing to join.
    """
    meshes = [m for m in meshes if m is not None and len(m.faces) > 0]
    if not meshes:
        return None

    all_v = []
    all_f = []
    for m in meshes:
        base = len(all_v)
        all_v.extend(np.asarray(m.vertices).tolist())
        for f in np.asarray(m.faces):
            all_f.append([int(f[0]) + base, int(f[1]) + base, int(f[2]) + base])

    # process=True merges coincident vertices, welding the shared cap/piece boundary
    joined = trimesh.Trimesh(np.array(all_v), np.array(all_f), process=True)
    # make winding consistent / normals outward, so a watertight result is a proper "volume"
    # (manifold3d booleans require this, and it fixes inverted/negative-volume pieces)
    joined.fix_normals()
    return joined


def spread_pieces(meshes, factor=0.5, center=None):
    """
    Move each mesh radially away from the common center by factor*(its center - overall center),
    for an exploded view (the Python twin of "Spread Pieces").  Returns new translated copies;
    the inputs are left untouched.
    """
    meshes = list(meshes)
    centers = [np.asarray(m.bounds).mean(axis=0) for m in meshes]
    if center is None:
        center = np.mean(centers, axis=0) if centers else np.zeros(3)
    else:
        center = np.asarray(center, dtype=float)[:3]

    out = []
    for m, c in zip(meshes, centers):
        moved = m.copy()
        moved.apply_translation(factor * (c - center))
        out.append(moved)
    return out


def mesh_boolean_fold(solid, features, signs=None):
    """
    Fold an ordered sequence of boolean operations onto a solid mesh -- the Python twin of the
    Grasshopper "Boolean Piece" component.

    solid: a trimesh.Trimesh (should be watertight; booleans on open meshes are unreliable).
    features: list of trimesh.Trimesh to boolean in, IN ORDER.  ("Feature" in the solid-modeling
              sense -- an ordered additive/subtractive operation on a body.)
    signs: list parallel to features; +1 = union, <=0 = subtract.  Defaults to all subtract.

    Order matters: each step acts on the result of the previous, e.g. signs [+1, -1, +1, -1]
    means union(f0), then subtract(f1), then union(f2), then subtract(f3).  Uses trimesh's
    exact 'manifold' backend (the manifold3d package), which is robust on clean manifolds.
    Returns the resulting trimesh.Trimesh.
    """
    try:
        import manifold3d  # noqa: F401  -- the exact boolean backend trimesh will use
    except ImportError as e:
        raise ImportError(
            "mesh booleans need the 'manifold3d' package (pip install manifold3d)") from e

    import warnings

    features = list(features)
    if signs is None:
        signs = [-1] * len(features)
    if len(signs) != len(features):
        raise ValueError("signs must be parallel to features")

    if not solid.is_watertight:
        warnings.warn("boolean solid is not watertight; the result may be wrong")

    result = solid.copy()
    for feature, sign in zip(features, signs):
        if sign > 0:
            result = trimesh.boolean.union([result, feature], engine='manifold')
        else:
            result = trimesh.boolean.difference([result, feature], engine='manifold')
    return result


class SurfacePiece():
    """
    A "Piece" of an algebraic surface.  Essentially, a union of Faces, with some additional interface.
    """

    def __init__(self, indices, surface):
        """ Initialize a SurfacePiece object with corresponding indices and surface

            :param indices: A list of nonsingular pieces' indices
            :param surface: Surface data
        """

        self.indices = indices
        self.surface = surface
        self.dimension = 2
        self.num_variables = surface.num_variables
        self.center = surface.center
        self.radius = surface.radius


        # memoized members:
        self.m_edge_pieces = None


    def __str__(self):
        """ toString method for SurfacePiece """
        result = "SurfacePiece with face indices:"
        result += "{}".format(self.indices)
        return result

    def __repr__(self):
        return str(self)

    def is_compact(self):
        """ Check whether a piece is:
            (1) compact (no edges touch the bounding sphere)
            (2) non-compact (at least 1 edge touches bounding sphere)


            Examples:
            sphere: (1 piece) - compact

            dingdong: (2 pieces) - one compact, one not compact

            octdong: (2 pieces) - both compact

            whitney: (2 pieces) - both non-compact
            paraboloid: (1 piece) - non compact

        """

        # bounding sphere
        sphere_curve = self.surface.sphere_curve.sampler_data

        for ii in self.indices:
            face = self.surface.faces[ii]
            if face['system top'] == "input_surf_sphere":
                return False

            if face['system bottom'] == "input_surf_sphere":
                return False

        # compact
        return True

    def centroid(self):
        """Compute the centroid of a piece"""

        def flatten_and_unique(list_nD):
            """helper fucntion for centroid to get a flat list of unique values"""
            return list(set([inner for outer in list_nD for inner in outer]))


        # the point indices on the piece: its samples when it has them (every face's, not
        # the last face's alone, as it was), else the raw faces' points
        unique_point_indices_this_piece = self.face_points(
            samples=self.surface.is_sampled(), as_indices=True, unique=True)
        #deref each point index to its point
        points = self.surface.extract_points()
        coordinates_this_piece = np.array([points[ind,:] for ind in unique_point_indices_this_piece])
        #return the mean as [x,y,z]
        return coordinates_this_piece.mean(axis=0)

    # point_singularities
    # the points on a piece ,  left and right edge will be degenerated
    # type critical

    def point_singularities(self):
        """ Compute the indices of the singularity points from a SurfacePiece object

            :rtype: A list of indices of point singularities
        """

        point_singularities = []

        for face_index in self.indices:
            surf = self.surface

            # if vertex type is singular, returns true
            face = surf.faces[face_index]

            curr_edge = -10

            # top - make this to function (later)
            # top is the index, system top is the where the index lives
            if(face['system top'] == 'input_critical_curve'):
                curr_edge = surf.critical_curve.edges[face['top']]
            elif(face['system top'] == 'input_surf_sphere'):
                curr_edge = surf.sphere_curve.edges[face['top']]
            else:
                for zz in range(len(surf.singular_curves)):
                    if(surf.singular_names[zz] == face['system top']):
                        curr_edge = surf.singular_curves[zz].edges[face['top']]

            # vertices
            for ii in range(3): # 0 is left, 1 is mid, 2 is right
                if(self.surface.vertices[curr_edge[ii]].is_of_type(VertexType.singular)):
                    point_singularities.append(curr_edge[ii])


            # check bottom is singular
            if(face['system bottom'] == 'input_critical_curve'):
                curr_edge = surf.critical_curve.edges[face['bottom']]
            elif(face['system bottom'] == 'input_surf_sphere'):
                curr_edge = surf.sphere_curve.edges[face['bottom']]
            else:
                for zz in range(len(surf.singular_curves)):
                    if(surf.singular_names[zz] == face['system bottom']):
                        curr_edge = surf.singular_curves[
                            zz].edges[face['bottom']]
            # vertices
            for ii in range(3):
                if(self.surface.vertices[curr_edge[ii]].is_of_type(VertexType.singular)):
                    point_singularities.append(curr_edge[ii])

            # now we check to the left
            for edge_ind in face['left']: # this thing itself is a list
                curr_edge = surf.critical_point_slices[face['middle slice index']].edges[edge_ind]
                for ii in range(3):
                    if(self.surface.vertices[curr_edge[ii]].is_of_type(VertexType.singular)):
                        point_singularities.append(curr_edge[ii])

                        # now we check to the right
            for edge_ind in face['right']: # this thing itself is a list
                curr_edge = surf.critical_point_slices[face['middle slice index']+1].edges[edge_ind]
                for ii in range(3):
                    if(self.surface.vertices[curr_edge[ii]].is_of_type(VertexType.singular)):
                        point_singularities.append(curr_edge[ii])

        return list(set(point_singularities))


    def write_skeleton_data(self):
        raise NotImplementedError

        pieces = self.separate_into_nonsingular_pieces()


        with open(self.generate_filename_no_ext("skeleton_piece")+".scad", "w") as f:
            pass


        return



    def plot(self, color, ax):
        return self.surface.plot(face_indices=self.indices, color=color, ax=ax)



    def _edges_touching(self):
        """
        computes dictionary of the curve edges on this piece of surface.  includes crit, sing, mid, and sphere.  
        this function will return a dict of sets of edge indices.  

        the intention of this function is to convert them to CurvePieces in a subsequent step.
        """

        from collections import defaultdict
        touching_curve_edge_indices = defaultdict(set)
    

        for face_index in self.indices:

            face = self.surface.faces[face_index]

            slice_ind = face['middle slice index']

            left = self.surface.critical_point_slices[slice_ind]
            for edge_ind in face['left']:
                touching_curve_edge_indices[left.inputfilename].add(edge_ind)
                

            right = self.surface.critical_point_slices[slice_ind+1]
            for edge_ind in face['right']:
                touching_curve_edge_indices[right.inputfilename].add(edge_ind)

            touching_curve_edge_indices[face['system top']].add(face['top'])
            touching_curve_edge_indices[face['system bottom']].add(face['bottom'])


            mid = self.surface.midpoint_slices[slice_ind]
            for ind, e in enumerate(mid.edges):
                if e[1]==face['midpoint']:
                    touching_curve_edge_indices[mid.inputfilename].add(ind)

        return dict(touching_curve_edge_indices) # change from a default_dict of sets to a dict of sets



    





    def edge_pieces(self):
        """
        takes a `dict` of `set`s of edge indices. 
        produces a `list` of `CurvePiece`s


        a kind of related note: the critical curve is very likely in the middle of the surface piece
        the boundary of a `SurfacePiece` is probably sphere or singular `CurvePiece`s.  
        it's possible the edges are degenerate, in case
        of nodal singularity.
        """

        
        # an act of memoization
        if self.m_edge_pieces:
            return self.m_edge_pieces


        touching_curve_edge_indices = self._edges_touching()

        self.m_edge_pieces = [] # will be a list of CurvePieces

        for curve_name, edge_indices in touching_curve_edge_indices.items():
            self.m_edge_pieces.extend( self.surface.curve_with_name(curve_name).break_into_pieces(edge_indices) )

        return self.m_edge_pieces



    def to_gh_dict(self, piece_index, include_smooth=True):
        """
        assemble this piece's data for the Grasshopper JSON export.

        meshes are expressed purely as triangle indices into the surface's unified vertex
        set; embedded curves as ordered vertex-index lists into the same set.  no vertex
        coordinates live here -- they are shared at the top level of the export.
        """

        mesh_raw = _mesh_triangles(self.surface.as_mesh_raw(self.indices))

        mesh_smooth = None
        if include_smooth and self.surface.is_sampled():
            # samples taken of the merged faces are asked for by the piece's merged faces, as in
            # export_smooth: as_mesh_smooth then reads its face indices as merged ones
            which = self.indices
            if getattr(self.surface, 'samples_of', 'faces') == 'merged':
                which = self.surface.merged_faces_of(self.indices)
            try:
                mesh_smooth = _mesh_triangles(self.surface.as_mesh_smooth(which))
            except br_except.SurfaceNotSampled:
                mesh_smooth = None

        curves = []
        for cp in self.edge_pieces():
            curves.append({
                "type": self.surface._curve_type_for_name(cp.curve.inputfilename),
                "curve_name": cp.curve.inputfilename,
                "vertex_indices": cp.to_point_indices(),
            })

        return {
            "piece_index": piece_index,
            "face_indices": list(self.indices),
            "mesh_smooth": mesh_smooth,
            "mesh_raw": mesh_raw,
            "curves": curves,
        }



    def generate_filename_no_ext(self,basename, ninds=3):
        """ 
        construct a filename for the piece, using face indices to make unique.  
        generates without an extension, so that it can be added later
        """

        return basename+'_'+ ('-'.join([str(i) for i in self.indices[: min(len(self.indices),ninds) ]])) # have to use `min` in case piece has <ninds faces on it


    def generate_filename_smooth(self, file_type=_default_file_type):
        return '{}.{}'.format( self.generate_filename_no_ext(basename=_default_piece_basename_smooth), file_type)

    def generate_filename_raw(self, file_type=_default_file_type):
        return '{}.{}'.format( self.generate_filename_no_ext(basename=_default_piece_basename_raw), file_type)


    def export_smooth(self, basename=_default_piece_basename_smooth,autoname_using_folder=False,file_type=_default_file_type):
        """Write the piece's smooth mesh; returns the file name.  Samples taken of the merged
        faces are asked for by the piece's merged faces (`Surface.merged_faces_of`), since
        `as_mesh_smooth` then reads its face indices as merged ones."""
        filename_no_ext = self.generate_filename_no_ext(basename)
        which = self.indices
        if getattr(self.surface, 'samples_of', 'faces') == 'merged':
            which = self.surface.merged_faces_of(self.indices)
        return self.surface.export_smooth(which,filename_no_ext,autoname_using_folder,file_type)


    def export_raw(self, basename=_default_piece_basename_raw,autoname_using_folder=False,file_type=_default_file_type):

        filename_no_ext = self.generate_filename_no_ext(basename)
        self.surface.export_raw(self.indices,filename_no_ext,autoname_using_folder,file_type)


    def as_mesh(self, smooth=None):
        """
        The `trimesh.Trimesh` for this piece.  smooth=None picks smooth when the surface is
        sampled, else raw (raw with raw, sampled with sampled).
        """
        if smooth is None:
            smooth = self.surface.is_sampled()
        if smooth:
            return self.surface.as_mesh_smooth(self.indices)
        return self.surface.as_mesh_raw(self.indices)


    def sphere_caps(self, smooth=None, resolution=4, tol=1e-3):
        """
        The faceted spherical cap mesh(es) closing this piece where it meets the bounding
        sphere.  See the module-level `sphere_cap_meshes`.  Returns a list of trimesh.
        """
        return sphere_cap_meshes(self.as_mesh(smooth), self.center, self.radius, resolution, tol)


    def flat_caps(self, smooth=None, resolution=1, tol=1e-3):
        """
        The faceted FLAT cap mesh(es) closing this piece where it meets the bounding sphere, with
        the apex at each loop's centroid.  See the module-level `flat_cap_meshes`.  Returns a list
        of trimesh.
        """
        return flat_cap_meshes(self.as_mesh(smooth), self.center, self.radius, resolution, tol)


    def as_closed_mesh(self, smooth=None, resolution=None, tol=1e-3, flat=False):
        """
        This piece joined with its cap(s) into a single welded (ideally watertight)
        `trimesh.Trimesh` -- the Rhino-free equivalent of (Sphere|Flat) Caps + Close Piece.  A
        piece bounded only by the sphere comes out watertight; one abutting a singular curve stays
        open there (check `.is_watertight`).

        flat=False uses spherical caps (hugging the sphere); flat=True uses flat fans to the loop
        centroid.  resolution defaults to 4 for spherical, 1 for flat.
        """
        if resolution is None:
            resolution = 1 if flat else 4
        mesh = self.as_mesh(smooth)
        if flat:
            caps = flat_cap_meshes(mesh, self.center, self.radius, resolution, tol)
        else:
            caps = sphere_cap_meshes(mesh, self.center, self.radius, resolution, tol)
        return join_meshes([mesh] + caps)


    def solidify_smooth(self, distance=_default_solidify_thickness, basename=_default_piece_basename_smooth, autoname_using_folder=False,file_type=_default_file_type):

        filename_no_ext = self.generate_filename_no_ext(basename)
        self.surface.solidify_smooth(distance, self.indices,filename_no_ext,autoname_using_folder,file_type)


    def solidify_raw(self, distance=_default_solidify_thickness, basename=_default_piece_basename_raw, autoname_using_folder=False,file_type=_default_file_type):

        filename_no_ext = self.generate_filename_no_ext(basename)
        self.surface.solidify_raw(distance, self.indices,filename_no_ext,autoname_using_folder,file_type)




    def face_points(self, samples=True, as_indices = False, unique = True):
        """
        Get the coordinates of the points on all the faces of the Piece of a Surface.

        if `samples`, then will return all samples on the Piece.  otherwise, will return the points of the raw faces.  

        - the computed point set should have no duplicates.  
        - i do not know what order the points will be in, sorry.
        """

        point_indices = list()

        if samples:
            # the samples may have been taken of the merged faces: then take this piece's
            for tri in self.surface.sample_triangles(
                    self.surface.merged_faces_of(self.indices)
                    if self.surface.samples_of == 'merged' else self.indices,
                    merged=self.surface.samples_of == 'merged'):
                point_indices.extend(tri)

        else:

            # need the midpoint of the face (comes from the mid of the mid), and the left/right/mid of all edges.
            touching_curve_edge_indices = self._edges_touching() # this is a dict of strings and lists-of-ints

            for curve_name, edge_indices in touching_curve_edge_indices.items():
                curve = self.surface.curve_with_name(curve_name)
                for e in edge_indices:
                    point_indices.extend(curve.edges[e])

        if unique:
            point_indices = list(set(point_indices))

        if as_indices:
            return point_indices

        else:# next, get the actual coordinates from the vertex set
            return self.surface.extract_points(indices=point_indices)



class Surface(Decomposition):
    """ Create a Surface object (Child class of Decomposition)

        :param Decomposition: Decomposition data from decomp file

    """

    def __init__(self, directory, is_embedded=False,embedded_into=None):
        """ Initialize a Surface Object

            :param directory: Directory of the surface folder
        """

        self.num_faces = 0
        self.num_midpoint_slices = 0
        self.num_critical_slices = 0
        self.num_singular_curves = 0
        self.singular_curve_multiplicities = []
        self.faces = {}    # stores all data from F.faces file
        self.midpoint_slices = []
        self.critical_point_slices = []
        self.critical_curve = []
        self.sphere_curve = []
        self.singular_curves = []
        self.singular_names = []
        # preimage (semicritical) curves: a bertini2-backend addition -- the
        # shadow arcs of the critical curve's image on the non-fold sheets.  Not
        # produced by classic bertini_real, so they are discovered by globbing
        # curve_preimage_* rather than counted in S.surf (backward compatible:
        # a classic decomposition simply has none).
        self.preimage_curves = []
        self.preimage_names = []
        # merged faces (F.merged, a bertini2-backend addition): None without them
        self.merged_faces = None
        self.edge_runs = {}
        self.sampler_data = []   # store all surface_sampler data
        self.samples_of = 'faces'   # what sampler_data is per: 'faces' or 'merged' faces

        Decomposition.__init__(self, directory, is_embedded,embedded_into)


        # automatically parse data files to gather curve data
        self.parse_surf(self.directory)
        self.gather_faces(self.directory)
        self.gather_merged_faces(self.directory)
        self.gather_curves(self.directory)
        try:
            self.gather_surface_samples(self.directory)
        except:
            print("no samples found")

    def __repr__(self):
        """ toString method for Surface """
        result = "Surface with:\n"
        result += f"{self.num_faces} faces\n"
        result += f"defined using {self.num_variables} variables\n"
        result += f"center, radius of sphere: {self.center}, {self.radius}\n"
        result += f"there are {self.num_critical_slices} crit slices\n"
        result += f"and {self.num_singular_curves} singular_curves with multiplicities {self.singular_curve_multiplicities}\n\n"
        result += f"computed point set has {len(self.vertices)} total points in it\n"
        result += f""
        return result

    def __str__(self):
        return repr(self)



    def is_sampled(self):
        """
        Query whether the surface has been sampled.  
        """

        return len(self.sampler_data) > 0


    # this should be a decorator
    def _require_samples(self):
        if not self.is_sampled():
            raise br_except.SurfaceNotSampled('trying to do something that requires samples, but surface is not sampled.  Sample and re-gather, or make do with raw / blocky data')



    def parse_surf(self, directory):
        """ Parse and store into surface data

            :param directory: Directory of the surface folder
        """
        surf_data = bertini_real.parse.parse_surf(directory)
        self.num_faces = surf_data[0]
        self.num_edges = surf_data[1]
        self.num_midpoint_slices = surf_data[2]
        self.num_critical_slices = surf_data[3]
        self.num_singular_curves = surf_data[4]
        self.singular_curve_multiplicities = surf_data[5]

    # def parse_vertex_types(self, directory):
    #     """ Parse and store vertex types data

    #     :param directory: Directory of the surface folder
    #     """
    #     vertex_types_data = bertini_real.parse.parse_vertex_types(directory)
    #     self.vertex_types_data = vertex_types_data

    def gather_faces(self, directory):
        """ Gather the faces of surface

            :param directory: Directory of the surface folder
        """
        self.faces = bertini_real.parse.parse_faces(directory)

    def gather_curves(self, directory):
        """ Gather the curves of surface

            :param directory: Directory of the surface folder
        """
        for ii in range(self.num_midpoint_slices):
            new_curve = Curve(directory + '/curve_midslice_' + str(ii),is_embedded=True,embedded_into=self)
            self.midpoint_slices.append(new_curve)
        for ii in range(self.num_critical_slices):
            new_curve = Curve(directory + '/curve_critslice_' + str(ii),is_embedded=True,embedded_into=self)
            self.critical_point_slices.append(new_curve)

        self.critical_curve = Curve(directory + '/curve_crit',is_embedded=True,embedded_into=self)
        self.sphere_curve = Curve(directory + '/curve_sphere',is_embedded=True,embedded_into=self)

        for ii in range(self.num_singular_curves):
            filename = directory + '/curve_singular_mult_' + \
                str(self.singular_curve_multiplicities[ii][0]) + '_' + str(
                    self.singular_curve_multiplicities[ii][1])
            new_curve = Curve(filename,is_embedded=True,embedded_into=self)
            self.singular_curves.append(new_curve)
            self.singular_names.append(new_curve.inputfilename)

        # preimage (semicritical) curves are not counted in S.surf (a bertini2
        # addition); discover them by globbing curve_preimage_N in index order
        import glob, os, re
        preimage_dirs = glob.glob(os.path.join(directory, 'curve_preimage_*'))
        def _preimage_index(path):
            m = re.search(r'curve_preimage_(\d+)$', path)
            return int(m.group(1)) if m else 0
        for d in sorted(preimage_dirs, key=_preimage_index):
            new_curve = Curve(d, is_embedded=True, embedded_into=self)
            self.preimage_curves.append(new_curve)
            self.preimage_names.append(new_curve.inputfilename)

    def gather_merged_faces(self, directory):
        """ Gather the merged faces of the surface, when the decomposition carries F.merged
            (a bertini2-backend addition, see `parse_merged_faces`): faces joined across the
            critical slices that cut nothing between them.  Each is a dict with the keys of a
            face in F.faces widened: 'top' and 'bottom' are lists of edge indices, one per
            face left to right; 'left' and 'right' are the sides on 'left slice index' and
            'right slice index'; plus 'faces' (its faces in F.faces), 'seams' and 'seeds'.
            `merged_faces` is None for a decomposition without them, and `edge_runs` holds
            per curve the runs of edges welded at the merged-away joints.

            :param directory: Directory of the surface folder
        """
        parsed = bertini_real.parse.parse_merged_faces(directory)
        if parsed is None:
            self.merged_faces, self.edge_runs = None, {}
            return
        merged, self.edge_runs = parsed
        for m in merged:
            faces = [self.faces[p] for p in m['faces']]
            for which in ('top', 'bottom'):
                systems = {f['system ' + which] for f in faces}
                if len(systems) != 1:
                    raise RuntimeError(f"merged face with faces {m['faces']} has its {which} "
                                       f"on several curves: {sorted(systems)}")
                m[which] = [f[which] for f in faces]
                m['system ' + which] = faces[0]['system ' + which]
            m['left'] = list(faces[0]['left'])
            m['right'] = list(faces[-1]['right'])
            m['num left'] = len(m['left'])
            m['num right'] = len(m['right'])
            m['seeds'] = [f['midpoint'] for f in faces]
        self.merged_faces = merged

    def _faces(self, merged):
        """The faces, or the merged faces when asked, which the decomposition must carry."""
        if not merged:
            return self.faces
        if self.merged_faces is None:
            raise RuntimeError('this decomposition carries no merged faces (no F.merged); '
                               'merge its bundle with python -m cellular.merge_faces and '
                               'export it with --merged')
        return self.merged_faces

    def merged_faces_of(self, face_indices):
        """The indices of the merged faces whose faces are among `face_indices` (a piece's
        faces, say).  A merged face lies within one piece, since no seam passes through a
        singular point; one with faces on both sides of the given set is an error."""
        wanted = set(face_indices)
        out = []
        for k, m in enumerate(self._faces(True)):
            inside = [p in wanted for p in m['faces']]
            if all(inside):
                out.append(k)
            elif any(inside):
                raise RuntimeError(f"merged face {k} has faces {m['faces']} only partly among "
                                   f"the given faces")
        return out

    def _collapsed_run_edges(self, curve_name, edge_indices):
        """The edges `edge_indices` of the named curve as (left, midpoint, right) triples, with
        each run of edges welded at merged-away joints (`edge_runs`) collapsed to one triple:
        the run's outer ends with the midpoint of its middle edge, the vertex the merge keeps."""
        edges = self.curve_with_name(curve_name).edges
        run_of = {}
        for run in self.edge_runs.get(curve_name, []):
            for k in run:
                run_of[k] = tuple(run)
        out = []
        group, group_run = [], None
        for k in list(edge_indices) + [None]:
            run = run_of.get(k) if k is not None else None
            if group and (k is None or run is None or run != group_run):
                if group_run is None:
                    out.extend([edges[j] for j in group])
                else:
                    keep = group_run if list(group_run) == group else tuple(group)
                    out.append([edges[group[0]][0], edges[keep[len(keep) // 2]][1],
                                edges[group[-1]][2]])
                group, group_run = [], None
            if k is None:
                break
            if k < 0:
                continue
            if run is None:
                out.append(edges[k])
            else:
                group.append(k)
                group_run = run
        return out

    def boundary_edges(self, face):
        """The edges around a face as (left, midpoint, right) vertex triples, each as the curve
        stores it, in the order the raw triangulation walks them: the top reversed, the bottom,
        the left side, the right side reversed; a degenerate edge is left out.  For a merged
        face (`merged_faces`) the top and bottom are their runs of face edges with the welded
        runs collapsed (`_collapsed_run_edges`), and the sides lie on its outer slices."""
        if isinstance(face['top'], list):
            top = self._collapsed_run_edges(face['system top'], face['top'])
            bottom = self._collapsed_run_edges(face['system bottom'], face['bottom'])
            left_slice, right_slice = face['left slice index'], face['right slice index']
        else:
            top = ([self.curve_with_name(face['system top']).edges[face['top']]]
                   if face['top'] >= 0 else [])
            bottom = ([self.curve_with_name(face['system bottom']).edges[face['bottom']]]
                      if face['bottom'] >= 0 else [])
            left_slice, right_slice = face['middle slice index'], face['middle slice index'] + 1
        left = [self.critical_point_slices[left_slice].edges[k] for k in face['left'] if k >= 0]
        right = [self.critical_point_slices[right_slice].edges[k] for k in face['right'] if k >= 0]
        around = ([list(reversed(e)) for e in top] + list(bottom) + list(left)
                  + [list(reversed(e)) for e in right])
        return [e for e in around if not (e[0] < 0 and e[1] < 0 and e[2] < 0)]

    def gather_surface_samples(self, directory):
        """ Gather the surface samples of surface: one list of triangles per face of F.faces,
            or, when the samples were taken of the merged faces (a bertini2-backend sampling of
            F.merged, one list per merged face), per merged face; `samples_of` says which,
            'faces' or 'merged'.

            :param directory: Directory of the surface folder
        """
        self.sampler_data = bertini_real.parse.parse_surface_samples(
            directory)
        self.samples_of = 'faces'
        if (self.merged_faces is not None and len(self.sampler_data) == len(self.merged_faces)
                and len(self.sampler_data) != self.num_faces):
            self.samples_of = 'merged'

    def sample_triangles(self, face_indices, merged=False):
        """The sampled triangles (global vertex index triples) of the given faces, or with
        `merged` of the given merged faces.  Samples taken of the faces serve a merged face
        as the union of its faces' samples; samples taken of the merged faces cannot be
        split into faces, so asking for faces then is an error."""
        self._require_samples()
        if merged == (self.samples_of == 'merged'):
            return [tri for ii in face_indices for tri in self.sampler_data[ii]]
        if merged:
            return [tri for ii in face_indices for p in self.merged_faces[ii]['faces']
                    for tri in self.sampler_data[p]]
        raise RuntimeError('the samples were taken of the merged faces, which do not split '
                           'into faces; ask for merged faces (merged=True)')



    def check_data(self):
        """ Check data """
        try:
            if self.dimension != 2:
                print('This function designed to work on surfaces decomposed with bertini_real.  your object has dimension ' + self.dimension)

        except:
            return

    def faces_nonsingularly_connected(self, seed_index):
        """ Compute the faces that are nonsingualrly connected

            :param seed_index: Index of seed
            :rtype: Two lists containing indices of connected and unconnected faces
        """
        self.check_data()

        new_indices = [seed_index]
        connected = []

        while not(new_indices == []):
            connected.extend(new_indices)
            new_indices = self.find_connected_faces(connected)

        connected.sort()
        set_num_faces = list(range(self.num_faces))

        unconnected = list(set(set_num_faces) - set(connected))

        return connected, unconnected

    def find_connected_faces(self, current):
        """ Find connected faces from current face

            :param current: Current face
            :rtype: List containing indices of connected faces

        """

        new_indices = []

        unexamined_indices = list(range(self.num_faces))

        unexamined_indices = list(set(unexamined_indices) - set(current))

        for ii in range(len(current)):
            c = current[ii]
            f = self.faces[c]  # unpack the current face
            deleteme = []

            for jj in range(len(unexamined_indices)):
                d = unexamined_indices[jj]
                g = self.faces[d]  # unpack the examined face

                if self.faces_nonsingularly_connect(f, g):
                    new_indices.append(d)
                    deleteme.append(d)

            unexamined_indices = list(set(unexamined_indices) - set(deleteme))

        return new_indices

    def faces_nonsingularly_connect(self, f, g):
        """ Check whether faces f and g are nonsingularly connected

            :param f: Current face
            :param g: Other face
            :rtype: Return True if f and g are nonsingularly connected, else False
        """
        val = False

        if self.cannot_possibly_meet(f, g):
            return val

        elif self.faces_meet_at_left(f, g):
            val = True

        elif self.faces_meet_at_right(f, g):
            val = True

        elif self.faces_meet_at_top(f, g):
            val = True

        elif self.faces_meet_at_bottom(f, g):
            val = True

        return val

    def cannot_possibly_meet(self, f, g):
        """ Check whether faces f and g cannot possibly meet (because they are in different fiber intervals of the projection)

            :param f: Current face
            :param g: Other face
            :rtype: Return True if f and g meet, else False
        """
        val = False

        if abs(f['middle slice index'] - g['middle slice index']) >= 2:
            val = True

        return val

    def faces_meet_at_left(self, f, g):
        """ Check whether faces f and g nonsingularly connected at left

            :param f: Current face
            :param g: Other face
            :rtype: Return True if f and g nonsingularly connected at left, else False
        """
        val = False

        for ii in range(f['num left']):
            e = self.critical_point_slices[
                f['middle slice index']].edges[f['left'][ii]]
            a = e[1]

            for jj in range(g['num left']):
                E = self.critical_point_slices[
                    g['middle slice index']].edges[g['left'][jj]]
                b = E[1]

                if a == b and not(is_edge_degenerate(e)) and not(is_edge_degenerate(E)):
                    val = True
                    return val

            for jj in range(g['num right']):
                E = self.critical_point_slices[
                    g['middle slice index'] + 1].edges[g['right'][jj]]
                b = E[1]

                if a == b and not(is_edge_degenerate(e)) and not(is_edge_degenerate(E)):
                    val = True
                    return val
        return val

    def faces_meet_at_right(self, f, g):
        """ Check whether faces f and g nonsingularly connected at right

            :param f: Current face
            :param g: Other face
            :rtype: Return True if f and g nonsingularly connected at right, else False
        """
        val = False

        for ii in range(f['num right']):
            e = self.critical_point_slices[
                f['middle slice index'] + 1].edges[f['right'][ii]]
            a = e[1]

            for jj in range(g['num left']):
                E = self.critical_point_slices[
                    g['middle slice index']].edges[g['left'][jj]]
                b = E[1]

                if a == b and not(is_edge_degenerate(e)) and not(is_edge_degenerate(E)):
                    val = True
                    return val

            for jj in range(g['num right']):
                E = self.critical_point_slices[
                    g['middle slice index'] + 1].edges[g['right'][jj]]
                b = E[1]

                if a == b and not(is_edge_degenerate(e)) and not(is_edge_degenerate(E)):
                    val = True
                    return val
        return val

    def faces_meet_at_top(self, f, g):
        """ Check whether faces f and g nonsingularly connected at top

            :param f: Current face
            :param g: Other face
            :rtype: Return True if f and g nonsingularly connected at top, else False
        """
        val = False

        if(f['system top'][0:15] == 'input_singcurve'):
            return val  # cannot meet singularly, because edge is singular
        else:
            # at least they are in the same interval
            if f['middle slice index'] != g['middle slice index']:
                return val

        if (f['system top'] == g['system top']):
            if (self.critical_curve.inputfilename == f['system top']):
                if (f['top'] == g['top']):
                    val = True
                    return val

        if (f['system top'] == g['system bottom']):
            if (self.critical_curve.inputfilename == f['system top']):
                if (f['top'] == g['bottom']):
                    val = True
                    return val

        return val

    def faces_meet_at_bottom(self, f, g):
        """ Check whether faces f and g nonsingularly connected at bottom

            :param f: Current face
            :param g: Other face
            :rtype: Return True if f and g nonsingularly connected at bottom, else False
        """
        val = False

        if(f['system bottom'][0:15] == 'input_singcurve'):
            return val  # cannot meet singularly, because edge is singular
        else:
            # at least they are in the same interval
            if f['middle slice index'] != g['middle slice index']:
                return val

        if (f['system bottom'] == g['system top']):
            if (self.critical_curve.inputfilename == f['system bottom']):
                if (f['bottom'] == g['top']):
                    val = True
                    return val

        if (f['system bottom'] == g['system bottom']):
            if (self.critical_curve.inputfilename == f['system bottom']):
                if (f['bottom'] == g['bottom']):
                    val = True
                    return val

        return val



    def separate_into_nonsingular_pieces(self):
        """
        Separate a surface into a list of pieces, connected at singularities.

        Two faces are in one piece when they share a nondegenerate edge of a critical slice (a
        side of each), or an edge of the critical curve in the same interval (the top or bottom
        of each), as `faces_nonsingularly_connect` decides for one pair; edges of the sphere and
        of singular curves join nothing.  Each face is indexed once by the edges it could share,
        and the pieces are the connected components, so the cost grows with the number of faces,
        not with its square: on a barth10 decomposition of 18636 faces, the pairwise search
        took over an hour.

        Pieces are ordered by their smallest face index, and each piece's face indices are
        sorted.
        """

        self.check_data()

        parent = list(range(self.num_faces))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        first_with = {}  # an edge a face could share -> the first face seen with it

        def join(face_index, shared):
            other = first_with.setdefault(shared, face_index)
            a, b = find(face_index), find(other)
            if a != b:
                parent[max(a, b)] = min(a, b)

        crit_name = self.critical_curve.inputfilename
        for ii, f in enumerate(self.faces):
            s = f['middle slice index']
            for slice_index, sides in ((s, f['left']), (s + 1, f['right'])):
                edges = self.critical_point_slices[slice_index].edges
                for k in sides or []:
                    e = edges[k]
                    if not is_edge_degenerate(e):
                        join(ii, ('slice edge midpoint', e[1]))
            for which in ('top', 'bottom'):
                if f['system ' + which] == crit_name:
                    join(ii, ('critical curve edge', s, f[which]))

        groups = defaultdict(list)
        for ii in range(self.num_faces):
            groups[find(ii)].append(ii)

        return [SurfacePiece(sorted(g), self) for g in sorted(groups.values(), key=min)]


    
    def curve_with_name(self, curve_name):

        if curve_name == self.critical_curve.inputfilename:
            return self.critical_curve

        if curve_name == self.sphere_curve.inputfilename:
            return self.sphere_curve

        for c in self.critical_point_slices:
            if curve_name == c.inputfilename:
                return c

        for c in self.midpoint_slices:
            if curve_name == c.inputfilename:
                return c 

        for c in self.singular_curves:
            if curve_name == c.inputfilename:
                return c

        for c in self.preimage_curves:
            if curve_name == c.inputfilename:
                return c

        raise RuntimeError(f'unable to find a curve with name {curve_name} in this surface')


    def _curve_type_for_name(self, curve_name):
        """
        classify an embedded curve by its `inputfilename` into one of the closed-vocabulary
        type tags used by the Grasshopper export.  match order mirrors `curve_with_name`.
        """

        if curve_name == self.critical_curve.inputfilename:
            return "critical"

        if curve_name == self.sphere_curve.inputfilename:
            return "sphere"

        for c in self.critical_point_slices:
            if curve_name == c.inputfilename:
                return "critslice"

        for c in self.midpoint_slices:
            if curve_name == c.inputfilename:
                return "midslice"

        if curve_name in self.singular_names:
            return "singular"

        return "unknown"


    def export_gh_json(self, filename="br_gh_export.json", include_smooth=True):
        """
        write a self-contained JSON describing this surface for the Grasshopper plugin.

        the file holds one unified vertex set (`vertices`); each nonsingular piece carries
        only triangle indices (raw and, when sampled, smooth) and the embedded curve pieces
        as ordered vertex-index lists -- all indices into the shared `vertices`.  this keeps
        the surface mesh and its embedded curves referring to the same points in Rhino.
        """

        # prime the extract_points memo cache with the full (no-arg) point set first, so
        # later per-piece mesh construction does not poison it with a partial set.
        points = self.extract_points()

        pieces = self.separate_into_nonsingular_pieces()

        contents = {
            "format_version": 2,
            "decomposition_type": "surface",
            "source_directory": self.directory,
            "num_variables": self.num_variables,
            "vertices": _points_to_xyz(points),
            "vertex_count": len(points),
            "sphere": self._sphere_dict(),
            "is_sampled": self.is_sampled(),
            "pieces": [p.to_gh_dict(ii, include_smooth) for ii, p in enumerate(pieces)],
        }

        # fold the singularity / connector data (locations, tangent-cone directions, parities)
        # into the same file, so Grasshopper has a single source instead of a second JSON.
        try:
            sing = self.singularity_connector_data()
        except Exception as e:
            print("WARNING: could not compute singularity connector data ({}); "
                  "exporting without singularities".format(e))
            sing = {"piece_names": [], "on_pieces": [], "locations": [],
                    "directions": [], "parities": []}
        contents["singularities"] = {
            "piece_names": sing["piece_names"],
            "locations": sing["locations"],
            "directions": sing["directions"],
            "parities": sing["parities"],
            "on_pieces": sing["on_pieces"],
        }

        # verify each piece's sphere curves are closed loops (they always should be, barring a
        # decomposition problem); warn loudly if not, so the issue is visible before Grasshopper.
        for pc in contents["pieces"]:
            for cv in pc["curves"]:
                if cv["type"] == "sphere":
                    vi = cv["vertex_indices"]
                    if len(vi) < 2 or vi[0] != vi[-1]:
                        print("WARNING: piece {} has a non-closed sphere curve ({}); "
                              "the decomposition may be incomplete".format(
                                  pc["piece_index"], cv["curve_name"]))

        with open(filename, "w") as f:
            json.dump(contents, f, indent=2)

        print("wrote " + filename)
        return filename


    def all_curves(self):

        the_curves = []

        the_curves.append(self.critical_curve)

        the_curves.append(self.sphere_curve)

        for c in self.singular_curves:
            the_curves.append(c)

        for c in self.critical_point_slices:
            the_curves.append(c)

        for c in self.midpoint_slices:
            the_curves.append(c)

        return the_curves


    def all_singular_points(self):
        """
        get absolutely all of the singular points.  
        """

        # there's a baked-in assumption that this surface is NOT contained in a higher-dimensional object.  this is valid right now because the top-dimensional thing Bertini_real can decompose is a surface.


        the_singularites = []
        for v in self.vertices:
            if v.is_of_type(VertexType.singular):
                the_singularites.append(v)

        return the_singularites

    def isolated_singularities(self):
        VertexType = bertini_real.vertex.VertexType

        the_singularites = []
        for v in self.vertices:
            if v.is_of_type(VertexType.singular) and v.is_of_type(VertexType.singular):
                the_singularites.append(v)

        return the_singularites


    def singularity_connector_data(self, pieces=None):
        """
        Compute the data needed to place plug/socket connectors at nodal singularities.

        For each nodal singularity that joins exactly two nonsingular pieces, find its
        location, the connector axis direction (the tangent-cone direction, from the Hessian
        of the defining polynomial at the singularity), and the per-piece parity (which side
        gets the plug vs the socket); also record which singularities lie on each piece.

        pieces: the pieces to consider, if already separated (all of them by default).  Only
        these count: a singularity where one of them meets a piece left out joins one piece
        here, not two, so it gets no connector.

        Needs the `bertini` parser and `sympy`, but only when there is at least one qualifying
        singularity.  Returns a dict of pure-Python (JSON-safe) values:
          { "piece_names": [str, ...],          # per piece
            "on_pieces":   [[int, ...], ...],    # per piece: compact singularity indices
            "locations":   [[x, y, z], ...],     # per singularity
            "directions":  [[x, y, z], ...],     # per singularity
            "parities":    [[int, ...], ...] }   # per singularity: a value per piece (-1/0/1)
        All lists are empty when there are no qualifying singularities.
        """

        pieces = self.separate_into_nonsingular_pieces() if pieces is None else list(pieces)
        piece_names = [p.generate_filename_smooth() for p in pieces]

        # nodal singularities and which pieces each is incident to
        sings_on_pieces = {}
        pieces_connected_to_sing = defaultdict(list)
        for ii, p in enumerate(pieces):
            sings_this_piece = p.point_singularities()  # indices into the vertex set
            sings_on_pieces[ii] = sings_this_piece
            for s in sings_this_piece:
                pieces_connected_to_sing[s].append(ii)

        # only singularities joining exactly two pieces receive a connector
        wanted = {k: v for k, v in pieces_connected_to_sing.items() if len(v) == 2}
        singindex2int = {s: i for i, s in enumerate(wanted.keys())}

        on_pieces = []
        for ii in range(len(pieces)):
            on_pieces.append([singindex2int[s] for s in wanted if s in sings_on_pieces[ii]])

        def unit_vector(vector):
            return vector / np.linalg.norm(vector)

        locations = []
        directions = []

        if wanted:
            import bertini as b2
            import sympy

            # the whole input file: b2's parser takes the CONFIG and INPUT sections as they are,
            # and the text after 'INPUT' alone ends in a bare END; that it rejects
            bsys = b2.parse.system(self.input)
            f = bsys.function(0)
            F = sympy.S(str(f).replace('unnamed_function', '').replace('function', '').replace('f', ''))
            variables = sorted(F.free_symbols, key=lambda s: s.name)
            H = sympy.hessian(F, variables)
            hessian_evalme = sympy.lambdify(variables, H, modules='numpy')

            for sing_index in wanted.keys():
                sing_coords = self.vertices[sing_index].point.real

                # tangent-cone direction: eigenvector of the Hessian belonging to the
                # odd-one-out (smallest) eigenvalue
                M = hessian_evalme(*sing_coords)
                q = np.linalg.eig(M)
                axis = np.real(q.eigenvectors[:, np.argmin(q.eigenvalues)])
                direction0 = unit_vector(np.asarray(axis, dtype=float))

                directions.append([float(x) for x in direction0])
                locations.append([float(x) for x in sing_coords])

        parities = [[0 for _ in range(len(pieces))] for _ in range(len(wanted))]
        for s, ps in wanted.items():
            parities[singindex2int[s]][ps[0]] = -1
            parities[singindex2int[s]][ps[1]] = 1

        return {
            "piece_names": piece_names,
            "on_pieces": on_pieces,
            "locations": locations,
            "directions": directions,
            "parities": parities,
        }


    def write_piece_data(self, pieces=None, compact_only=False):
        """
        Writes the piece data for the plugs and sockets (br_surf_piece_data.scad and .json, in
        the working directory): the pieces' smooth mesh file names, the nodes at which two of the
        pieces meet, with a direction and a location for each, and which side of each node each
        piece is on.

        pieces: the pieces to write, if already separated (all of them by default).
        compact_only: only the pieces that do not touch the bounding sphere.  A node where a kept
            piece meets a left-out one is then no longer a node between two pieces, so it gets no
            plug and socket: the finite pieces can be printed and joined on their own.

        Returns the list of pieces written, in the order the files index them.
        """

        pieces = self.separate_into_nonsingular_pieces() if pieces is None else list(pieces)
        if compact_only:
            pieces = [p for p in pieces if p.is_compact()]
        data = self.singularity_connector_data(pieces)
        piece_names = data["piece_names"]
        singularities_on_pieces = data["on_pieces"]
        sing_directions_as_list = data["directions"]
        sing_locations_as_list = data["locations"]
        parity_of_sing_by_piece = data["parities"]
        allPoints = []

        #open and auto write the data(piece file names (without extensions), all sings of pieces, sing directions in order of sing index, sing coords in order of sing index) of the piece
        with open("br_surf_piece_data.scad", "w") as f:
            f.write(f'piece_names = [')
            f.write('"{}"'.format( '","'.join(piece_names) ))

            f.write('];\n')
            f.write(f'singularities_on_pieces = {singularities_on_pieces};\n')
            f.write(f'sing_directions = {sing_directions_as_list};\n')
            f.write(f'sing_locations = {sing_locations_as_list};\n')

            f.write(f'parities = {parity_of_sing_by_piece};\n')
            f.write(f'conn_size = 0.01;\n') #hard coded, but needs to be automatically computed
        print('br_surf_piece_data.scad')

        # Option 2: custom JSON encoder
        class NumpyEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                if isinstance(obj, np.integer):
                    return int(obj)
                if isinstance(obj, np.floating):
                    return float(obj)
                return super().default(obj)



        #open and auto write piece data to a json file
        with open("br_surf_piece_data.json", "w") as j:
            j.write(json.dumps({"piece_names": piece_names,
            "singularities_on_pieces": singularities_on_pieces,
            "sing_directions": sing_directions_as_list,
            "sing_locations": sing_locations_as_list,
            "parities" : parity_of_sing_by_piece},indent=4,cls=NumpyEncoder))
        print('wrote br_surf_piece_data.json')


        # with open("centroids.json", "w") as c:
        #     for centroid in centroids:
        #         c.write(str(centroid)+"\n")
        # print('wrote centroids.json')



        with open("allPoints.json", "w") as a:
            for point in allPoints:
                a.write("\n".join([str(s) for s in point]) + "\n")
        print('wrote allPoints.json')
        return pieces


    def as_mesh_smooth(self, which_faces=None, keep_all_vertices=True, merged=None):
        """
        Compute a `Trimesh` object from the `trimesh` library for the corresponding faces using sampled data.  Raises if the surface is not sampled.

        which_faces: either None for all faces, or a list-like of ints indicating the indices of the surface faces you want.
        keep_all_vertices: bool, by default True.  Unused vertices will be kept or merged.  This value influences `Trimesh`'s `process` parameter.
        merged: `which_faces` indexes the merged faces (`merged_faces`) rather than the faces; by default, whichever the samples were taken of (`samples_of`).  See `sample_triangles`.

        See https://trimsh.org/trimesh.html#trimesh.Trimesh.
        """
        if not self.sampler_data:
            raise br_except.SurfaceNotSampled('no surface samples found.  run sampler, or re-gather and pickle')
        if merged is None:
            merged = self.samples_of == 'merged'
        if which_faces is None:
            which_faces = range(len(self._faces(merged)))

        points = self.extract_points()
        triangles = [[tri[0], tri[1], tri[2]]
                     for tri in self.sample_triangles(which_faces, merged=merged)]

        should_trimesh_process = False if keep_all_vertices==True else True

        A = trimesh.Trimesh(np.array(points), np.array(triangles), process=should_trimesh_process)
        A.fix_normals()
        return A


    def as_mesh_raw(self, which_faces=None, keep_all_vertices=True, merged=False):
        """
        Compute a `Trimesh` object from the `trimesh` library for the corresponding faces using raw (unsmoothed or blocky) data: every edge around a face makes two triangles with the face's midpoint (`boundary_edges`).

        which_faces: either None for all faces, or a list-like of ints indicating the indices of the surface faces you want.
        keep_all_vertices: bool, by default True.  Unused vertices will be kept or merged.  This value influences `Trimesh`'s `process` parameter.
        merged: triangulate the merged faces (`merged_faces`, a decomposition carrying F.merged) instead of the faces; `which_faces` then indexes them.  A merged face spanning many slices becomes a few big triangles about its seed, the facets of a cut stone, where its faces were many slivers.

        See https://trimsh.org/trimesh.html#trimesh.Trimesh.
        """
        faces = self._faces(merged)
        if which_faces is None:
            which_faces = range(len(faces))

        points = self.extract_points()

        triangles = []
        for ii in which_faces:
            face = faces[ii]
            if not merged and face['middle slice index'] == -1:
                continue
            # fan each edge around the face to the face midpoint; skip degenerate triangles (a
            # repeated vertex index, e.g. from a degenerate curve edge).  these are zero-area,
            # and if kept they make the mesh non-manifold -- an (a, a, mid) face contributes the
            # {a, mid} edge twice.
            for e in self.boundary_edges(face):
                for tri in ((e[0], e[1], face['midpoint']), (e[1], e[2], face['midpoint'])):
                    if len(set(tri)) == 3:
                        triangles.append(list(tri))

        # honor keep_all_vertices like as_mesh_smooth: process=False keeps the full global
        # vertex set so the faces index into extract_points() (the unified set), and does not
        # merge coincident-but-distinct vertices (which would fuse sheets at singularities).
        should_trimesh_process = False if keep_all_vertices == True else True
        raw_mesh = trimesh.Trimesh(np.array(points), np.array(triangles),
                                   process=should_trimesh_process)
        raw_mesh.fix_normals()
        return raw_mesh








    def export_pieces_glb(self, filename='br_pieces.glb', coloring='contrast', compact_only=False, flip=False, pieces=None, saturation=0.9, value=0.8, markers=True, piece_colors=None, facets='faces', group=None):
        """
        Write the surface's nonsingular pieces to one GLB file (glTF), each piece its own named mesh
        (`piece_017_compact`, `piece_142_cut`, ...) with its own colored material, for Blender and
        other viewers, and a .json beside it naming every piece's faces, compactness, centroid,
        color and nodes.

        Every piece points outward with uniform normals, exact from the input's polynomial when the
        surface is one function of three variables; see `bertini_real.surface.pieces_glb`.

        coloring: 'contrast' (every piece its own hue, pieces meeting at a node far apart in hue),
            'azimuth' (hue by the angle of a piece's centroid about the z axis), 'height' (by its
            z) or 'spread' (by piece order).
        compact_only: only the pieces that do not touch the bounding sphere.
        flip: point every piece inward instead.
        pieces: the pieces to write, if already separated.
        saturation, value: the colors' HSV saturation and value, as seen (sRGB); the defaults are
            deep jewel tones, and a lower saturation gives pastels.
        markers: a separate mesh named `nodes`, a small gold sphere at every node where two
            pieces meet.
        piece_colors: one (r, g, b) in [0, 1] per piece, as seen (sRGB), in the order of
            `pieces`; overrides coloring, saturation and value.
        facets: 'faces' triangulates every face of F.faces; 'merged' triangulates the merged
            faces (F.merged) instead, so each stone has a few big facets rather than slivers;
            'smooth' uses the surface's samples (samp.surfsamp), when it has them.
        group: the name of the root node every piece is parented under, and of the glTF scene
            (by default the file's name): in Blender, one empty, or with 4.2 and later one
            collection, whose eye shows or hides the whole subject.

        Returns a summary: pieces, triangles, how many point outward, the fraction of the area that
        faces its normals, and whether the normals are exact.
        """
        from bertini_real.surface import pieces_glb
        return pieces_glb.export(self, filename, coloring=coloring, compact_only=compact_only,
                                 flip=flip, pieces=pieces, saturation=saturation, value=value,
                                 markers=markers, piece_colors=piece_colors, facets=facets,
                                 group=group)


    def export_raw(self, which_faces=None, basename=_default_surface_basename_raw, autoname_using_folder=False, file_type=_default_file_type, keep_all_vertices=True):
        """ 
        Export raw decomposition of surface

        returns the name of the file which was saved
        """

        mesh = self.as_mesh_raw(which_faces,keep_all_vertices)
        return export_mesh(mesh, basename, autoname_using_folder, file_type)


    def export_smooth(self, which_faces=None, basename=_default_surface_basename_smooth, autoname_using_folder=False, file_type=_default_file_type, keep_all_vertices=True):
        """ 
        Export smooth decomposition of surface
        
        returns the name of the file which was saved
        """

        mesh = self.as_mesh_smooth(which_faces,keep_all_vertices)
        return export_mesh(mesh, basename, autoname_using_folder, file_type)









    def solidify_raw(self, distance=_default_solidify_thickness, which_faces=None, basename=_default_surface_basename_raw+'solidified', autoname_using_folder=False, file_type=_default_file_type, keep_all_vertices=True):
        """
        Solidify raw version of surface.

        Available formats include {'stl', 'obj'}.
        Default file format given by `_default_file_type`

        returns the name of the file which was saved
        """

        mesh = self.as_mesh_raw(which_faces,keep_all_vertices)
        solid = solidify_mesh(mesh,distance)
        return export_mesh(solid, basename, autoname_using_folder, file_type)




    def solidify_smooth(self, distance=_default_solidify_thickness, which_faces=None, basename=_default_surface_basename_raw+'solidified', autoname_using_folder=False, file_type=_default_file_type, keep_all_vertices=True):
        """
        Solidify smooth version of surface.  Requires that the surface has been sampled using `sampler`

        Available formats include {'stl', 'obj'}.
        Default file format given by `_default_file_type`

        returns the name of the file which was saved
        """

        mesh = self.as_mesh_smooth(which_faces,keep_all_vertices)
        solid = solidify_mesh(mesh,distance)
        return export_mesh(solid, basename, autoname_using_folder, file_type)







    #################


    #  _______  _        ______     _______  _______    _______           _______  _______  _______  _______  _______
    # (  ____ \( (    /|(  __  \   (  ___  )(  ____ \  (  ____ \|\     /|(  ____ )(  ____ \(  ___  )(  ____ \(  ____ \
    # | (    \/|  \  ( || (  \  )  | (   ) || (    \/  | (    \/| )   ( || (    )|| (    \/| (   ) || (    \/| (    \/
    # | (__    |   \ | || |   ) |  | |   | || (__      | (_____ | |   | || (____)|| (__    | (___) || |      | (__
    # |  __)   | (\ \) || |   | |  | |   | ||  __)     (_____  )| |   | ||     __)|  __)   |  ___  || |      |  __)
    # | (      | | \   || |   ) |  | |   | || (              ) || |   | || (\ (   | (      | (   ) || |      | (
    # | (____/\| )  \  || (__/  )  | (___) || )        /\____) || (___) || ) \ \__| )      | )   ( || (____/\| (____/\
    # (_______/|/    )_)(______/   (_______)|/         \_______)(_______)|/   \__/|/       |/     \|(_______/(_______/



    #############
