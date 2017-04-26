import ParsingFunctions
import os

class Decomposition(object):
    def parse_decomp(self, directory):
		decomposition_data = ParsingFunctions.parse_decomposition(directory)
		self.inputfilename = decomposition_data['input file name']
		self.pi = decomposition_data['Pi info']
		self.patch = decomposition_data['Patch Vectors']
		self.radius = decomposition_data["radius"]
		self.center = decomposition_data["center"]
		self.num_patches = decomposition_data['num patches']

    def read_input(self, directory):
		filename = directory + '/' + self.inputfilename
		if not os.path.isfile(filename):
			print "Could not find input file in current directory: %s" %(directory)
		else:
			with open(filename, 'r') as f:
				self.input = f.read()
