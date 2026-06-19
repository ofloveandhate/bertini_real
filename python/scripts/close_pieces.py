# Rhino-free pipeline: split a surface decomposition into its nonsingular pieces, cap each
# where it meets the bounding sphere, join into (ideally watertight) solids, optionally spread
# them apart for viewing, and export each as an STL.  This mirrors the Grasshopper components
# (Sphere Caps + Close Piece + Spread Pieces) for people who don't have a Rhino license.
#
# Watertight solids require a *sampled* decomposition: the raw (blocky) mesh is non-manifold,
# so raw pieces close visually but are not watertight.  The script reports this honestly.
#
# Usage (from a decomposition's working directory, or pass the output_dim_2_comp_0 folder):
#   python close_pieces.py [folder] [--resolution N] [--spread F] [--smooth | --raw]
#                          [--basename NAME] [--combined]

import argparse
import sys

import bertini_real as br
from bertini_real.surface import Surface, spread_pieces, join_meshes


def main(argv=None):
    ap = argparse.ArgumentParser(description="Close surface pieces with sphere caps and export STLs (no Rhino needed).")
    ap.add_argument("folder", nargs="?", default=None,
                    help="path to an output_dim_2_comp_0 folder; omit to gather() from the current directory")
    ap.add_argument("--resolution", type=int, default=4, help="radial subdivisions of each sphere cap (default 4)")
    ap.add_argument("--spread", type=float, default=0.0, help="exploded-view factor; 0 = pieces in place (default 0)")
    ap.add_argument("--smooth", dest="smooth", action="store_true", default=None, help="force sampled (smooth) meshes")
    ap.add_argument("--raw", dest="smooth", action="store_false", help="force raw (blocky) meshes")
    ap.add_argument("--basename", default="br_closed_piece", help="output STL basename (default br_closed_piece)")
    ap.add_argument("--combined", action="store_true", help="also export one combined STL of all pieces")
    args = ap.parse_args(argv)

    if args.folder:
        surface = Surface(args.folder)
    else:
        decomposition = br.data.gather()
        if not isinstance(decomposition, Surface):
            sys.exit("this script is for surface (dimension 2) decompositions")
        surface = decomposition

    pieces = surface.separate_into_nonsingular_pieces()
    print("{} nonsingular piece(s)".format(len(pieces)))

    closed = []
    for i, piece in enumerate(pieces):
        mesh = piece.as_closed_mesh(smooth=args.smooth, resolution=args.resolution)
        closed.append(mesh)
        if mesh is None:
            print("  piece {}: empty".format(i))
        else:
            print("  piece {}: {} faces, watertight={}".format(i, len(mesh.faces), mesh.is_watertight))

    # report the singularity / connector data too (locations + tangent-cone directions)
    try:
        sing = surface.singularity_connector_data()
        print("{} nodal singularity connector(s)".format(len(sing["locations"])))
    except Exception as e:
        print("(singularity connector data unavailable: {})".format(e))

    meshes = [m for m in closed if m is not None]
    if args.spread:
        meshes = spread_pieces(meshes, factor=args.spread)
        print("spread pieces apart by factor {}".format(args.spread))

    for i, mesh in enumerate(meshes):
        outname = "{}_{}.stl".format(args.basename, i)
        mesh.export(outname)
        print("wrote {}".format(outname))

    if args.combined and meshes:
        combined = join_meshes(meshes)
        if combined is not None:
            combined.export(args.basename + "_all.stl")
            print("wrote {}_all.stl".format(args.basename))


if __name__ == "__main__":
    main()
