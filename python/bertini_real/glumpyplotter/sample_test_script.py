#!/usr/local/bin/python3
import bertini_real as br


def function(x, y, z):
    return x + y + z


fn = function

br.glumpyplotter.plot_surface_samples()
# br.glumpyplotter.plot_surface_samples(cmap='inferno')
# br.glumpyplotter.plot_surface_samples(cmap='inferno', color_function=fn)
# br.glumpyplotter.plot_critical_curve(cmap='inferno', color_function=fn)
