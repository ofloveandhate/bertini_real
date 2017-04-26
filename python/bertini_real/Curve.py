import ParsingFunctions
from Decomposition import Decomposition

class Curve(Decomposition):
	def __init__(self, directory):
		self.num_edges = 0
		self.edges = []
		self.directory = directory
		self.input = None
		self.inputfilename = None
		self.num_variables = 0
		self.dimension = 0
		self.pi = []
		self.num_patches = 0
		self.patch  = {}
		self.radius = 0
		self.center_size = 0
		self.center = []
		self.sampler_data = []

		# automatically parse data files to gather curve data
		self.parse_decomp(self.directory)
		self.parse_edge(self.directory)
		self.parse_curve_sampler(self.directory)
		self.read_input(self.directory)


	def parse_edge(self, directory):
		edge_data = ParsingFunctions.parse_Edges(directory)
		self.num_edges = edge_data['number of edges']
		self.edges = edge_data['edges']

	def parse_curve_sampler(self, directory):
		self.sampler_data = ParsingFunctions.parse_Curve_Sampler(directory)
