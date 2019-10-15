
import bertini_real.data
import bertini_real.surface
import bertini_real.curve
import bertini_real.util
import bertini_real.plot
import bertini_real.glumpyplotter
import bertini_real.tmesh
import bertini_real.anaglypy


def gather_and_plot():
    bertini_real.data.gather()
    return bertini_real.plot.plot()


class vertex_types():

    def __init__(self):
        self.unset = 0
        self.critical = 1
        self.semicritical = 2
        self.midpoint = 4
        self.isolated = 8
        self.new = 16
        self.curve_sample_point = 32
        self.surface_sample_point = 64
        self.problematic = 256

    def unset(self):
        return self.__unset

    def critical(self):
        return self.__critical

    def semicritical(self):
        return self.__semicritical

    def midpoint(self):
        return self.__midpoint

    def isolated(self):
        return self.__isolated

    def new(self):
        return self.__new

    def curve_sample_point(self):
        return self.__curve_sample_point

    def surface_sample_point(self):
        return self.__surface_sample_point

    def problematic(self):
        return self.__problematic
