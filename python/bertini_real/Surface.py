import ParsingFunctions
import os
class Surface(object):
	def __init__(self, directory):
		self.directory = directory
		self.inputfilename = None
		self.num_variables = 0
		self.dimension = 9
		self.pi = []
		self.num_patches = 0
		self.patch = {}
		self.radius = 0
		self.center_size = 0
		self.center = []
		self.num_faces = 0
		self.num_edges = 0
		self.num_midpoint_slices = 0
		self.num_critical_slices = 0
		self.num_singular_curves = 0
		self.singular_curve_multiplicities = []
		self.faces = {}    # stores all data from F.faces file
		self.midpoint_slices = []
		self.critical_point_slices = []
		self.critical_curve = []
		self.sphere_curve = []
		self.singular_curves = []
		self.singular_names = []
		self.surface_sampler_data = []   # store all surface_sampler data
		self.input = None

		# automatically parse data files to gather curve data
		self.parse_decomp(self.directory)
		self.parse_surf(self.directory)
		self.gather_faces(self.directory)
		self.gather_curves(self.directory)
		self.gather_surface_samples(self.directory)
		self.read_input(self.directory)

	def parse_decomp(self, directory):
		decomposition_data = ParsingFunctions.parse_decomposition(directory)
		self.inputfilename = decomposition_data['input file name']
		self.pi = decomposition_data['Pi info'] 
		self.patch = decomposition_data['Patch Vectors']
		self.radius = decomposition_data["radius"]
		self.center = decomposition_data["center"]
		self.num_patches = decomposition_data['num patches']
	def parse_surf(self, directory):
		surf_data = ParsingFunctions.parse_Surf(directory)
		self.num_faces = surf_data[0]
		self.num_edges = surf_data[1]
		self.num_midpoint_slices = surf_data[2]
		self.num_critical_point_slices = surf_data[3]
		self.num_singular_curves = surf_data[4]
		self.singular_curve_multiplicities = surf_data[5]

	def gather_faces(self, directory):
		self.faces = ParsingFunctions.parse_Faces(directory)
	def gather_curves(self, directory):
		for ii in xrange(self.num_midpoint_slices):
			new_curve = Curve(directory + '/curve_midslice_' + str(ii))
			self.midpoint_slices.append(new_curve)
		for ii in xrange(self.num_critical_slices):
			new_curve = Curve(directory + '/curve_critslice_' + str(ii))
			self.critical_point_slices.append(new_curve)
		
		critical_curve = Curve(directory + '/curve_crit')
		self.critical_curve.append(critical_curve)
		sphere_curve = Curve(directory + '/curve_crit')

		for ii in xrange(self.num_singular_curves):
			filename = directory + '/curve_singular_mult_' + str(self.singular_curve_multiplicities[ii][0]) +'_' + str(self.singular_curve_multiplicities[ii][1])
			new_curve = Curve(filename)
			self.singular_curves.append(new_curve)
			self.singular_names.append(new_curve.inputfilename)


	def gather_surface_samples(self, directory):
		self.surface_sampler_data = ParsingFunctions.parse_Surface_Sampler(directory)

	def read_input(self, directory):
		filename = directory + '/' + self.inputfilename
		if not os.path.isfile(filename):
			print "Could not find input file in current directory: %s" %(directory)
		else:
			with open(filename, 'r') as f:
				self.input = f.read()









class Curve(Surface):
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
		self.curve_sampler_data = ParsingFunctions.parse_Curve_Sampler(directory)
	







