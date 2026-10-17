
# inspired by https://github.com/pypa/warehouse/
from bertini_real.__about__ import (
    __summary__,
    __title__,
    __uri__,
    __version__,
    __version_info__,
    __author__,
    __email__
)

    # __commit__,
    # __copyright__,
    # __license__,



import bertini_real.data

import bertini_real.edge
import bertini_real.face

import bertini_real.curve
import bertini_real.surface

import bertini_real.util
import bertini_real.plot

import bertini_real.paths

# bertini_real.glumpyplotter is deprecated (2026-10): it needs glumpy and a working OpenGL
# stack, which has become too hard to get running.  It is no longer imported here; import it
# explicitly if you have glumpy installed, and expect it to be removed.  Use
# bertini_real.plot (matplotlib) or Surface.export_pieces_glb for a viewer such as Blender.

import bertini_real.anaglypy


def gather_and_plot():
    decomposition = bertini_real.data.gather()
    return bertini_real.plot.plot(decomposition)
