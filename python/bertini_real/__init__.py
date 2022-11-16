
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
import bertini_real.surface
import bertini_real.curve
import bertini_real.util
import bertini_real.plot
try:
    import bertini_real.glumpyplotter
except ImportError as e:
    print(f'unable to import bertini_real.glumpyplotter.  if you want to use opengl-accelerated rendering of surfaces, install the module `glumpy` (and probably an opengl library, too)')

import bertini_real.anaglypy


def gather_and_plot():
    decomposition = bertini_real.data.gather()
    return bertini_real.plot.plot(decomposition)
