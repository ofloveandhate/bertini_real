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
import os
import re
import struct

import numpy as np
import trimesh

#: a vertex whose |grad F| is below this fraction of the median has no usable normal (a node)
NODE_GRADIENT_FRACTION = 1e-6

#: the colorings :func:`hues` knows
COLORINGS = ('contrast', 'azimuth', 'height', 'spread')

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


def neighbours_at_nodes(node_lists):
    """{piece: set of pieces} joined by a shared node, from each piece's list of nodes."""
    at = {}
    for k, nodes in enumerate(node_lists):
        for n in nodes:
            at.setdefault(n, []).append(k)
    out = {k: set() for k in range(len(node_lists))}
    for ks in at.values():
        for a in ks:
            out[a].update(b for b in ks if b != a)
    return out


def contrasting_hues(neighbours, count, seed=0, sweeps=400):
    """One hue per piece, every piece its own (``count`` hues evenly around the circle), placed so
    pieces that meet at a node are far apart in hue: simulated annealing over the assignment,
    swapping two pieces' hues, against a soft minimum of the hue distances across shared nodes.
    Deterministic for a given ``seed``."""
    import random
    rng = random.Random(seed)
    hue = [k / count for k in range(count)]
    rng.shuffle(hue)
    sigma = 0.08                                       # about 30 degrees

    def distance(a, b):
        d = abs(hue[a] - hue[b]) % 1.0
        return min(d, 1.0 - d)

    def penalty(k):
        return sum(math.exp(-distance(k, j) / sigma) for j in neighbours[k])

    steps = sweeps * count
    temperature, cooling = 1.0, (1e-4 / 1.0) ** (1.0 / steps)   # from 1 down to 1e-4
    for _ in range(steps):
        a, b = rng.randrange(count), rng.randrange(count)
        temperature *= cooling
        if a == b:
            continue
        before = penalty(a) + penalty(b)
        hue[a], hue[b] = hue[b], hue[a]
        after = penalty(a) + penalty(b)
        if after > before and rng.random() >= math.exp((before - after) / temperature):
            hue[a], hue[b] = hue[b], hue[a]
    return hue


def closest_neighbour_hues(hue, neighbours):
    """The smallest hue distance, in degrees, between any two pieces that meet at a node."""
    worst = 180.0
    for a, bs in neighbours.items():
        for b in bs:
            d = abs(hue[a] - hue[b]) % 1.0
            worst = min(worst, 360.0 * min(d, 1.0 - d))
    return worst


def hues(centroids, coloring, neighbours=None):
    """One hue in [0, 1) per piece: ``contrast`` every piece its own hue, pieces meeting at a node
    far apart (:func:`contrasting_hues`), ``azimuth`` by the angle of its centroid about the z axis
    (a rainbow sweeping around), ``height`` by its centroid's z (red at the bottom to violet at the
    top), ``spread`` by piece order, each 0.618 of a turn from the last (neighbours in the list
    contrast)."""
    if coloring == 'contrast':
        return contrasting_hues(neighbours, len(centroids))
    if coloring == 'azimuth':
        return [math.atan2(c[1], c[0]) / (2 * math.pi) % 1.0 for c in centroids]
    if coloring == 'height':
        z = np.array([c[2] for c in centroids])
        return list((z - z.min()) / max(z.max() - z.min(), 1e-300) * 0.85)
    if coloring == 'spread':
        return [(k * _GOLDEN) % 1.0 for k in range(len(centroids))]
    raise ValueError(f"coloring must be one of {COLORINGS}, not {coloring!r}")


def srgb_to_linear(channel):
    """One sRGB channel in [0, 1] as linear light: glTF's baseColorFactor is linear, and a color
    chosen in sRGB and written as it is comes out washed toward white."""
    if channel <= 0.04045:
        return channel / 12.92
    return ((channel + 0.055) / 1.055) ** 2.4


def signed_volume(mesh, about):
    """The signed volume swept from ``about`` to the mesh: for a closed, consistently wound mesh,
    its enclosed volume, positive exactly when the normals point out."""
    v, f = np.asarray(mesh.vertices) - np.asarray(about, dtype=float), np.asarray(mesh.faces)
    a, b, c = (v[f[:, j]] for j in range(3))
    return float(np.einsum('ij,ij->i', np.cross(b - a, c - a), a).sum()) / 6.0


def oriented_piece_mesh(surface, piece, gradient=None, flip=False, facets='panels'):
    """A piece's raw mesh, its own vertices only, wound consistently and pointing outward, with
    exact normals when ``gradient`` is given (see the module docstring).  ``facets='merged'``
    triangulates the piece's merged faces instead of its panels."""
    import contextlib
    import io
    merged = facets == 'merged'
    which = surface.merged_faces_of(piece.indices) if merged else piece.indices
    with contextlib.redirect_stdout(io.StringIO()):
        mesh = surface.as_mesh_raw(which_faces=which, keep_all_vertices=False, merged=merged)
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


def node_markers(points, radius):
    """One mesh of small spheres, one at each of ``points``."""
    sphere = trimesh.creation.icosphere(subdivisions=2, radius=radius)
    return trimesh.util.concatenate([sphere.copy().apply_translation(p) for p in points])


def name_glb_scene(filename, name):
    """Give the one scene in a written GLB the name ``name``, in place.  A GLB is a 12-byte
    header, a JSON chunk and a binary chunk; the JSON chunk is rewritten (padded to four bytes
    with spaces, as the format asks) and the header's total length updated.  Blender 4.2 and
    later import a glTF scene as a collection of that name."""
    with open(filename, 'rb') as f:
        data = f.read()
    magic, version, _length = struct.unpack_from('<III', data, 0)
    if magic != 0x46546C67:
        raise ValueError(f"{filename} is not a GLB")
    json_length, json_type = struct.unpack_from('<II', data, 12)
    if json_type != 0x4E4F534A:
        raise ValueError(f"{filename}: the first chunk is not JSON")
    tree = json.loads(data[20:20 + json_length].decode('utf-8'))
    for scene in tree.get('scenes', []):
        scene['name'] = name
    text = json.dumps(tree, separators=(',', ':')).encode('utf-8')
    text += b' ' * (-len(text) % 4)
    rest = data[20 + json_length:]
    body = struct.pack('<II', len(text), json_type) + text + rest
    with open(filename, 'wb') as f:
        f.write(struct.pack('<III', magic, version, 12 + len(body)) + body)


def export(surface, filename, coloring='contrast', compact_only=False, flip=False, pieces=None,
           saturation=0.9, value=0.8, roughness=0.55, markers=True, marker_radius=None,
           piece_colors=None, facets='panels', group=None):
    """Write the surface's nonsingular pieces to ``filename`` (a .glb) and a .json beside it naming
    every piece: its faces, whether it is compact, its centroid, its color, and the nodes on it.
    Colors are chosen in HSV with ``saturation`` and ``value`` (deep jewel tones by default; lower
    the saturation for pastels) as sRGB, and written as the linear colors glTF expects, so a viewer
    shows the color chosen.  The JSON gives each color as sRGB hex.  With ``markers``, a separate
    mesh named ``nodes`` holds a small gold sphere at every node where two pieces meet
    (``marker_radius``, by default 0.6% of the bounding sphere's radius), to toggle in a viewer.
    ``piece_colors``, one (r, g, b) in [0, 1] per piece as seen (sRGB), in the order of
    ``pieces``, overrides ``coloring``, ``saturation`` and ``value``.  ``facets='merged'``
    triangulates each piece's merged faces (``Surface.merged_faces``) instead of its panels.
    Everything is parented under one root node named ``group`` (by default the file's name
    without its extension), which Blender imports as an empty whose eye icon, Shift-clicked
    in the outliner, shows or hides the whole subject; the glTF scene gets the same name, so
    Blender 4.2 and later import the subject as a collection of that name, whose eye does the
    same.  Several subjects imported into one scene can then share its lights and camera.
    Returns a summary: pieces, triangles, how many point outward, the area fraction that faces its
    normals, whether the normals are exact, and for ``contrast`` the closest hues across a node.
    See the module docstring."""
    if facets not in ('panels', 'merged'):
        raise ValueError(f"facets must be 'panels' or 'merged', not {facets!r}")
    if group is None:
        group = re.sub(r'\.glb$', '', os.path.basename(filename))
    if pieces is None:
        pieces = surface.separate_into_nonsingular_pieces()
    gradient = gradient_from_input(surface.input)
    meshes = [oriented_piece_mesh(surface, p, gradient, flip, facets) for p in pieces]
    centroids = [np.asarray(m.centroid) for m in meshes]
    node_lists = [sorted(int(n) for n in p.point_singularities()) for p in pieces]
    neighbours = neighbours_at_nodes(node_lists)
    if piece_colors is not None:
        if len(piece_colors) != len(pieces):
            raise ValueError(f"piece_colors has {len(piece_colors)} colors for {len(pieces)} pieces")
        coloring = 'given'
        chosen = [tuple(float(t) for t in c[:3]) for c in piece_colors]
        piece_hues = [colorsys.rgb_to_hsv(*c)[0] for c in chosen]
    else:
        piece_hues = hues(centroids, coloring, neighbours)
        chosen = [colorsys.hsv_to_rgb(h, saturation, value) for h in piece_hues]

    scene = trimesh.Scene(base_frame=group)
    kept, record = [], []
    for k, (piece, mesh, centroid, srgb) in enumerate(zip(pieces, meshes, centroids, chosen)):
        compact = piece.is_compact()
        if compact_only and not compact:
            continue
        rgba = [int(round(255 * srgb_to_linear(t))) for t in srgb] + [255]
        hex_color = '#%02x%02x%02x' % tuple(int(round(255 * t)) for t in srgb)
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
                           color=hex_color,
                           nodes=node_lists[k], triangles=int(len(mesh.faces))))

    shared = sorted({n for k, ns in enumerate(node_lists) for n in ns
                     if not compact_only or pieces[k].is_compact()
                     if sum(n in other for other in node_lists) >= 2})
    if markers and shared:
        radius = marker_radius or 0.006 * float(surface.radius)
        points = [np.real(np.asarray(surface.vertices[n].point[:3])).astype(float) for n in shared]
        gold = trimesh.visual.material.PBRMaterial(
            name='nodes', baseColorFactor=[255, 196, 64, 255], metallicFactor=1.0,
            roughnessFactor=0.25)
        dots = node_markers(points, radius)
        dots = trimesh.Trimesh(dots.vertices, dots.faces, vertex_normals=dots.vertex_normals,
                               visual=trimesh.visual.TextureVisuals(material=gold), process=False)
        scene.add_geometry(dots, node_name='nodes', geom_name='nodes')

    scene.export(filename, include_normals=True)
    name_glb_scene(filename, group)
    json_name = re.sub(r'\.glb$', '', filename) + '.json'
    with open(json_name, 'w') as f:
        json.dump(dict(decomposition=surface.directory, coloring=coloring,
                       compact_only=compact_only, flipped=flip, facets=facets, group=group,
                       exact_normals=gradient is not None, pieces=record), f, indent=1)
    summary = check([m for m, _, _ in kept], [c for _, _, c in kept],
                    np.asarray(surface.center, dtype=float))
    summary.update(triangles=sum(r['triangles'] for r in record), facets=facets,
                   exact_normals=gradient is not None, filename=filename, json=json_name,
                   node_markers=len(shared) if markers else 0,
                   closest_hues_across_a_node=closest_neighbour_hues(piece_hues, neighbours))
    return summary
