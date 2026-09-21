"""
Regenerate the figures for tutorials/capping_and_joining.rst from the dingdong decomposition.

Run from the repo with the package importable, e.g.:
    PYTHONPATH=python python python/docs/tutorials/capping_and_joining_pictures/make_images.py \
        --dingdong test/surface/dingdong/output_dim_2_comp_0
Needs: trimesh, manifold3d (mesh booleans), bertini (tangent-cone directions), matplotlib.
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import trimesh
from bertini_real.surface import Surface, mesh_boolean_fold, spread_pieces

HERE = os.path.dirname(os.path.abspath(__file__))


def render(items, outpath, elev=22, azim=-60):
    """items: list of (trimesh, rgb-tuple)."""
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111, projection="3d")
    allpts = []
    for mesh, color in items:
        V = np.asarray(mesh.vertices)
        F = np.asarray(mesh.faces)
        tris = V[F]
        pc = Poly3DCollection(tris)
        shade = 0.45 + 0.55 * np.clip(mesh.face_normals[:, 2], 0, 1)
        cols = np.clip(np.array(color)[None, :] * shade[:, None], 0, 1)
        pc.set_facecolor(cols)
        pc.set_edgecolor((0, 0, 0, 0.10))
        pc.set_linewidth(0.1)
        ax.add_collection3d(pc)
        allpts.append(V)
    P = np.vstack(allpts)
    mins, maxs = P.min(0), P.max(0)
    c, r = (mins + maxs) / 2.0, (maxs - mins).max() / 2.0
    ax.set_xlim(c[0] - r, c[0] + r)
    ax.set_ylim(c[1] - r, c[1] + r)
    ax.set_zlim(c[2] - r, c[2] + r)
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, outpath), dpi=130)
    plt.close(fig)
    print("wrote", outpath)


def square_prism(width, length, loc, direction):
    box = trimesh.creation.box(extents=(width, width, length))  # centered at origin, along Z
    box.apply_transform(trimesh.geometry.align_vectors([0, 0, 1.0], direction))
    box.apply_translation(loc)
    return box


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dingdong", default="test/surface/dingdong/output_dim_2_comp_0")
    args = ap.parse_args()

    s = Surface(args.dingdong)
    pieces = s.separate_into_nonsingular_pieces()
    blue = (0.36, 0.56, 0.85)
    gold = (0.85, 0.66, 0.30)

    # --- flat vs spherical caps, on the piece that meets the sphere ---
    # the big "cone" piece has a large base cap; a low side view shows the dome-vs-flat silhouette
    cap_piece = max(range(len(pieces)), key=lambda i: pieces[i].as_mesh().bounding_box.volume)
    render([(pieces[cap_piece].as_closed_mesh(flat=False), blue)], "dingdong_sphere_caps.png", elev=8)
    render([(pieces[cap_piece].as_closed_mesh(flat=True), blue)], "dingdong_flat_caps.png", elev=8)

    # --- square hole in one piece, matching (smaller) rod in the other ---
    sd = s.singularity_connector_data()
    loc = np.array(sd["locations"][0], dtype=float)
    direction = np.array(sd["directions"][0], dtype=float)
    direction = direction / np.linalg.norm(direction)
    parity = sd["parities"][0]
    socket_idx = parity.index(-1)   # gets the hole
    plug_idx = parity.index(1)      # gets the rod

    # size the rod/hole to the SMALL (socket) piece so the hole leaves a ring of material
    width = 0.30
    length = 1.4                   # centered on the singularity; reaches into both pieces
    clearance = 0.85               # rod a bit smaller than the hole, so it slides in

    socket = pieces[socket_idx].as_closed_mesh()
    plug = pieces[plug_idx].as_closed_mesh()

    hole = square_prism(width, length, loc, direction)
    rod = square_prism(width * clearance, length, loc, direction)

    socket_holed = mesh_boolean_fold(socket, [hole], [-1])
    plug_rodded = mesh_boolean_fold(plug, [rod], [+1])

    render([(socket_holed, blue)], "dingdong_square_hole.png")
    render([(plug_rodded, gold)], "dingdong_square_rod.png")

    # assembled (rod seated in the hole) and exploded
    render([(socket_holed, blue), (plug_rodded, gold)], "dingdong_joined.png")
    exploded = spread_pieces([socket_holed, plug_rodded], factor=0.8)
    render([(exploded[0], blue), (exploded[1], gold)], "dingdong_exploded.png")

    print("watertight: socket_holed={}, plug_rodded={}".format(
        socket_holed.is_watertight, plug_rodded.is_watertight))


if __name__ == "__main__":
    main()
