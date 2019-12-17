"""
.. module:: decomposition
    :platform: Unix, Windows
    :synopsis: The Decomposition contains methods to read input file and parse decomposition.

"""
import bertini_real.parse
import os

class Decomposition(object):
    """ Create a Decomposition object (Parent Class) based on the decomposition

        :param object: Curve/Surface object

    """

    def __init__(self, directory, is_embedded=False):
        """ Initialize a Decomposition object

            :param directory: directory path
            :param is_embedded
        """
        self.directory = directory
        self.is_embedded = is_embedded
        self.input = None
        self.inputfilename = None
        self.num_variables = 0
        self.pi = []
        self.num_patches = 0
        self.patch = {} 
        self.radius = 0 
        self.center_size = 0 
        self.center = []
        self.parse_decomp(self.directory)
        self.read_input(self.directory)
        if not self.is_embedded:
            self.vertices, self.filenames = bertini_real.data.gather_vertices(self.directory)

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
