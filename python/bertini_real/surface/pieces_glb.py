"""
    :platform: Unix, Windows
    :synopsis: A surface's nonsingular pieces as one GLB file (glTF), each piece its own named mesh
               with its own material, with uniform outward normals, for Blender and other viewers.

Use it through :meth:`bertini_real.surface.Surface.export_pieces_glb`.

**Orientation.**  Every piece points OUTWARD, so a viewer can use the normals as written.  A
compact piece is closed (pinched only at its nodes), and its signed volume is positive exactly
when its normals point out of the region it bounds; a piece that comes out negative is turned over
as a whole.  A piece cut by the bounding sphere is not closed, and outward means away from the
sphere's center: the same signed sum, about the center.  ``flip=True`` points every piece inward.

**Normals.**  For a surface given by one polynomial F in three variables, read from the
decomposition's own input file, a vertex's normal is +-grad F / |grad F|, exact, with ONE sign per
piece: grad F / |grad F| is a single continuous normal field over a piece (connected, with no node
inside it), so it agrees with the outward orientation everywhere on the piece or nowhere.  Every
triangle is wound to face that field; the raw cells meet at T-junctions, so the mesh's own edges
cannot orient a piece.  So a faceted raw mesh shades smooth.  At a node grad F = 0 and there is no
normal; a vertex there takes the area-weighted normal of its own piece's triangles around it, so
each piece meeting at the node gets its own.  Any other surface gets the mesh's normals, oriented
through its edges as well as they allow.
"""

import colorsys
import json
import math
import re

import numpy as np
import trimesh

#: a vertex whose |grad F| is below this fraction of the median has no usable normal (a node)
NODE_GRADIENT_FRACTION = 1e-6

#: the colorings :func:`hues` knows
COLORINGS = ('azimuth', 'height', 'spread')

_GOLDEN = (math.sqrt(5) - 1) / 2


def gradient_from_input(input_text):
    """The gradient of a surface given as ONE function of THREE variables in a Bertini input file,
    as a numpy function of (x, y, z) arrays returning an (n, 3) array; None when the input is not
    one function of three variables, or does not parse.

    Reads the INPUT section: variable groups, constants and subfunctions (assignments to names that
    are not functions, substituted in order), and the function.  ``^`` is power, ``I`` the
    imaginary unit and ``Pi`` pi, as in Bertini."""
    import sympy
    from sympy.parsing.sympy_parser import parse_expr

    if input_text is None:
        return None
    match = re.search(r'\bINPUT\b(.*?)\bEND\s*;', input_text, re.S)
    if not match:
        return None
    statements = [s.strip() for s in match.group(1).split(';') if s.strip()]
    variables, functions, values = [], [], {}
    local = {'I': sympy.I, 'Pi': sympy.pi}
    try:
        for s in statements:
            words = s.split(None, 1)
            if words[0] in ('variable_group', 'hom_variable_group'):
                variables += [v.strip() for v in words[1].split(',')]
                local.update({v: sympy.Symbol(v) for v in variables})
            elif words[0] == 'function':
                functions += [f.strip() for f in words[1].split(',')]
            elif words[0] in ('constant', 'parameter', 'pathvariable', 'random', 'random_real'):
                continue
            elif '=' in s:
                name, expression = (t.strip() for t in s.split('=', 1))
                value = parse_expr(expression.replace('^', '**'),
                                   local_dict=dict(local, **values))
                values[name] = value
    except Exception:
        return None
    if len(variables) != 3 or len(functions) != 1 or functions[0] not in values:
        return None
    symbols = [local[v] for v in variables]
    F = values[functions[0]]
    partials = sympy.lambdify(symbols, [sympy.diff(F, s) for s in symbols], 'numpy')

    def gradient(points):
        points = np.asarray(points, dtype=float)
        columns = partials(points[:, 0], points[:, 1], points[:, 2])
        return np.column_stack([np.broadcast_to(np.real(c), (len(points),)) for c in columns])
    return gradient


def hues(centroids, coloring):
    """One hue in [0, 1) per piece: ``azimuth`` by the angle of its centroid about the z axis (a
    rainbow sweeping around), ``height`` by its centroid's z (red at the bottom to violet at the
    top), ``spread`` by piece order, each 0.618 of a turn from the last (neighbours in the list
    contrast)."""
    if coloring == 'azimuth':
        return [math.atan2(c[1], c[0]) / (2 * math.pi) % 1.0 for c in centroids]
    if coloring == 'height':
        z = np.array([c[2] for c in centroids])
        return list((z - z.min()) / max(z.max() - z.min(), 1e-300) * 0.85)
    if coloring == 'spread':
        return [(k * _GOLDEN) % 1.0 for k in range(len(centroids))]
    raise ValueError(f"coloring must be one of {COLORINGS}, not {coloring!r}")


def signed_volume(mesh, about):
    """The signed volume swept from ``about`` to the mesh: for a closed, consistently wound mesh,
    its enclosed volume, positive exactly when the normals point out."""
    v, f = np.asarray(mesh.vertices) - np.asarray(about, dtype=float), np.asarray(mesh.faces)
    a, b, c = (v[f[:, j]] for j in range(3))
    return float(np.einsum('ij,ij->i', np.cross(b - a, c - a), a).sum()) / 6.0


def oriented_piece_mesh(surface, piece, gradient=None, flip=False):
    """A piece's raw mesh, its own vertices only, wound consistently and pointing outward, with
    exact normals when ``gradient`` is given (see the module docstring)."""
    import contextlib
    import io
    with contextlib.redirect_stdout(io.StringIO()):
        mesh = surface.as_mesh_raw(which_faces=piece.indices, keep_all_vertices=False)
    mesh.remove_unreferenced_vertices()
    about = mesh.centroid if piece.is_compact() else np.asarray(surface.center, dtype=float)
    if gradient is None:
        trimesh.repair.fix_winding(mesh)
        if (signed_volume(mesh, about) < 0) != flip:
            mesh.invert()
        return trimesh.Trimesh(mesh.vertices, mesh.faces, vertex_normals=mesh.vertex_normals,
                               process=False)

    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces).copy()
    g = gradient(vertices)
    size = np.linalg.norm(g, axis=1)
    usable = size > NODE_GRADIENT_FRACTION * float(np.median(size))
    normals = np.zeros_like(g)
    normals[usable] = g[usable] / size[usable, None]

    a, b, c = (vertices[faces[:, j]] for j in range(3))
    backwards = np.einsum('ij,ij->i', np.cross(b - a, c - a), normals[faces].sum(axis=1)) < 0
    faces[backwards] = faces[backwards][:, [0, 2, 1]]
    oriented = trimesh.Trimesh(vertices, faces, process=False)
    normals[~usable] = np.asarray(oriented.vertex_normals)[~usable]

    if (signed_volume(oriented, about) < 0) != flip:
        faces = faces[:, [0, 2, 1]]
        normals = -normals
    return trimesh.Trimesh(vertices, faces, vertex_normals=normals, process=False)


def check(meshes, compact, center):
    """For meshes about to be written: how many point outward (signed volume about their centroid
    if compact, about ``center`` if cut), and the fraction of their area that faces its normals."""
    outward = 0
    facing_area = total_area = 0.0
    for mesh, is_compact in zip(meshes, compact):
        about = mesh.centroid if is_compact else center
        outward += signed_volume(mesh, about) > 0
        n = np.asarray(mesh.vertex_normals)
        v, f = np.asarray(mesh.vertices), np.asarray(mesh.faces)
        a, b, c = (v[f[:, j]] for j in range(3))
        cross = np.cross(b - a, c - a)
        facing = np.einsum('ij,ij->i', cross, n[f].sum(axis=1)) > 0
        area = 0.5 * np.linalg.norm(cross, axis=1)
        facing_area += float(area[facing].sum())
        total_area += float(area.sum())
    return dict(pieces=len(meshes), outward=int(outward),
                area_facing_normals=facing_area / total_area if total_area else 1.0)


def export(surface, filename, coloring='azimuth', compact_only=False, flip=False, pieces=None,
           saturation=0.78, value=0.95, roughness=0.55):
    """Write the surface's nonsingular pieces to ``filename`` (a .glb) and a .json beside it naming
    every piece: its faces, whether it is compact, its centroid, its color, and the nodes on it.
    Returns a summary: pieces, triangles, how many point outward, the area fraction that faces its
    normals, and whether the normals are exact.  See the module docstring."""
    if pieces is None:
        pieces = surface.separate_into_nonsingular_pieces()
    gradient = gradient_from_input(surface.input)
    meshes = [oriented_piece_mesh(surface, p, gradient, flip) for p in pieces]
    centroids = [np.asarray(m.centroid) for m in meshes]
    piece_hues = hues(centroids, coloring)

    scene = trimesh.Scene()
    kept, record = [], []
    for k, (piece, mesh, centroid, hue) in enumerate(zip(pieces, meshes, centroids, piece_hues)):
        compact = piece.is_compact()
        if compact_only and not compact:
            continue
        r, g, b = colorsys.hsv_to_rgb(hue, saturation, value)
        rgba = [int(round(255 * t)) for t in (r, g, b)] + [255]
        name = f"piece_{k:03d}" + ('_compact' if compact else '_cut')
        material = trimesh.visual.material.PBRMaterial(
            name=name, baseColorFactor=rgba, metallicFactor=0.0, roughnessFactor=roughness,
            doubleSided=True)
        colored = trimesh.Trimesh(mesh.vertices, mesh.faces, vertex_normals=mesh.vertex_normals,
                                  visual=trimesh.visual.TextureVisuals(material=material),
                                  process=False)
        scene.add_geometry(colored, node_name=name, geom_name=name)
        kept.append((colored, name, compact))
        record.append(dict(piece=k, name=name, compact=compact, faces=list(piece.indices),
                           centroid=[float(c) for c in centroid],
                           color='#%02x%02x%02x' % tuple(rgba[:3]),
                           nodes=sorted(int(n) for n in piece.point_singularities()),
                           triangles=int(len(mesh.faces))))

    scene.export(filename, include_normals=True)
    json_name = re.sub(r'\.glb$', '', filename) + '.json'
    with open(json_name, 'w') as f:
        json.dump(dict(decomposition=surface.directory, coloring=coloring,
                       compact_only=compact_only, flipped=flip,
                       exact_normals=gradient is not None, pieces=record), f, indent=1)
    summary = check([m for m, _, _ in kept], [c for _, _, c in kept],
                    np.asarray(surface.center, dtype=float))
    summary.update(triangles=sum(r['triangles'] for r in record),
                   exact_normals=gradient is not None, filename=filename, json=json_name)
    return summary
