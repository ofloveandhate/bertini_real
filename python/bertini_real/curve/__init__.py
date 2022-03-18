"""
    :platform: Unix, Windows
    :synopsis: The Curve is a child class of Decomposition. Curve contains methods that parse curve edge and samples data.
"""
import bertini_real.parse
from bertini_real.decomposition import Decomposition


class Curve(Decomposition):
    """ Create a Curve object (Child class of Decomposition)

        :param Decomposition: Decomposition data from decomp file

    """
    def __init__(self, directory, is_embedded=False,embedded_into=None):
        """ Initialize a Curve Object

            :param directory: Directory of the curve folder
        """
        
        self.num_edges = 0
        self.edges = []
        self.sampler_data = None

        Decomposition.__init__(self, directory, is_embedded,embedded_into)


        # automatically parse data files to gather curve data
        self.parse_edge(self.directory)
        try:
            self.parse_curve_samples(self.directory)
        except FileNotFoundError:
            print("no samples to gather") 

    def parse_edge(self, directory):
        """ Parse and store curve edges data

            :param directory: Directory of the curve folder

        """
        edge_data = bertini_real.parse.parse_edges(directory)
        self.num_edges = edge_data['number of edges']
        self.edges = edge_data['edges']

    def parse_curve_samples(self, directory):
        """ Parse and store curve samples data

            :param directory: Directory of the curve folder

        """
        self.sampler_data = bertini_real.parse.parse_curve_samples(directory)

    def __str__(self):
        """ toString method for Curve """
        result = "curve with:\n"
        result += "{} edges".format(self.num_edges)
        return result
