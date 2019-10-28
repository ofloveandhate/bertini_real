
import bertini_real.data
import bertini_real.surface
import bertini_real.curve
import bertini_real.util
import bertini_real.plot
import bertini_real.glumpyplotter
import bertini_real.tmesh
import bertini_real.anaglypy
import bertini_real.vertextype


def gather_and_plot():
    bertini_real.data.gather()
    return bertini_real.plot.plot()
