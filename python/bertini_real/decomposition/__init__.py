"""
    :platform: Unix, Windows
    :synopsis: The Decomposition contains methods to read input file and parse decomposition.
"""
import bertini_real.parse
import numpy as np
import os

class Decomposition(object):
    """ Create a Decomposition object (Parent Class) based on the decomposition

        :param object: Curve/Surface object

    """

    def __init__(self, directory, is_embedded=False, embedded_into=None):
        """ Initialize a Decomposition object

            :param directory: directory path
            :param is_embedded
            :param 
        """
        self.directory = directory

        self.is_embedded = is_embedded
        self.embedded_into = embedded_into

        self.input = None
        self.inputfilename = None
        self.num_variables = 0
        self.pi = []
        self.num_patches = 0
        self.patch = {} 
        self.radius = 0 
        self.center_size = 0 
        self.center = []
        self.dimension = 0

        self.br_version = bertini_real.__version_info__

        self.parse_decomp(self.directory)
        self.read_input(self.directory)

        self._memoized_data = {}

        if not self.is_embedded:
            self.vertices, self.filenames = bertini_real.data.gather_vertices(self.directory)
        else:
            if self.embedded_into is None:
                raise br_except.EmbeddedIssue("parameter `embedded_into` cannot be unset if the decomposition is embedded")

            self.vertices = self.embedded_into.vertices
            

    def parse_decomp(self, directory):
        """ Parse and store decomposition data

            :param directory: Directory of the decomposition

        """
        decomposition_data = bertini_real.parse.parse_decomposition(directory)
        self.inputfilename = decomposition_data['input file name']
        self.pi = decomposition_data['pi info']
        self.patch = decomposition_data['patch vectors']
        self.radius = decomposition_data["radius"]
        self.center = decomposition_data["center"]
        self.num_patches = decomposition_data['num patches']
        self.num_variables = decomposition_data['num_variables']
        self.dimension = decomposition_data['dimension']

    # parse edge here

    def read_input(self, directory):
        """ Read input file

            :param directory: Directory of the input

        """
        filename = directory + '/' + self.inputfilename
        if not os.path.isfile(filename):
            print("Could not find input file in current directory: %s" % (directory))
        else:
            with open(filename, 'r') as f:
                self.input = f.read()



    def extract_points(self, indices=None):
        """ Helper method
            Extract points from vertices as a list

            :rtype: numpy 2d array.  

        """

        if indices is None:
            indices = np.arange(len(self.vertices))

        if self.is_embedded:
            return self.embedded_into.extract_points()

        if '_memoized_data' in dir(self) and 'points' in self._memoized_data:
            return self._memoized_data['points']



        points = []

        for ii in indices:
            vertex = self.vertices[ii] # unpack via a reference

            point = [None] * self.num_variables

            for jj in range(self.num_variables):
                point[jj] = vertex.point[jj].real
            points.append(point)

        points = np.array(points)

        self._memoized_data['points'] = points.real

        return points
