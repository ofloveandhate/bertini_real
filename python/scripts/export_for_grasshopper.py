# Generate the intermediary JSON consumed by the bertini_real Grasshopper components.
#
# This writes one self-contained file (default: br_gh_export.json) holding a single unified
# vertex set; for a surface, each nonsingular piece is a mesh (triangle indices) plus the
# curves embedded on it (vertex-index polylines); for a curve, each piece is a polyline.
# Point the "Surface Read GH JSON" (or "Curve Read GH JSON") Grasshopper component at it.
#
# Usage:
#   # from a decomposition's working directory (the folder containing `Dir_Name`):
#   python export_for_grasshopper.py [output.json]
#
#   # or point it explicitly at an output_dim_X_comp_Y folder from anywhere:
#   python export_for_grasshopper.py path/to/output_dim_2_comp_0 [output.json]
#
# Works for both curves (dim 1) and surfaces (dim 2).

import os
import re
import sys

import bertini_real as br
from bertini_real.curve import Curve
from bertini_real.surface import Surface


def export_from_directory(directory, filename):
    """Build the decomposition directly from an output_dim_X_comp_Y folder and export it."""
    basename = os.path.basename(os.path.normpath(directory))
    match = re.search(r"dim_(\d+)", basename)
    if not match:
        raise SystemExit(
            f"could not determine the dimension from folder name {basename!r}; "
            "expected something like output_dim_2_comp_0"
        )

    dimension = int(match.group(1))
    if dimension == 1:
        decomposition = Curve(directory)
    elif dimension == 2:
        decomposition = Surface(directory)
    else:
        raise SystemExit(f"dimension {dimension} not supported (only curves=1 and surfaces=2)")

    return decomposition.export_gh_json(filename)


def main(argv):
    directory = None
    filename = "br_gh_export.json"

    # a positional arg that is an existing directory is the decomposition folder;
    # any other positional arg is the output filename.
    for arg in argv[1:]:
        if os.path.isdir(arg):
            directory = arg
        else:
            filename = arg

    if directory is not None:
        out = export_from_directory(directory, filename)
    else:
        # idiomatic path: gather from the current working directory (needs `Dir_Name`)
        out = br.data.gather_and_export_gh(filename)

    # export_gh_json already prints the filename it wrote; point at the next step.
    print(f"open {os.path.abspath(out)} with the 'Surface Read GH JSON' / 'Curve Read GH JSON' component")


if __name__ == "__main__":
    main(sys.argv)
