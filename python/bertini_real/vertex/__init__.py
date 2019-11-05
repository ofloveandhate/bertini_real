"""
.. module:: vertex
    :platform: Unix, Windows
    :synopsis: This module contains different types of VertexType

"""
# create getters for the Vertex class
# change all the vertex lookup to the new vertex syntax
# remove all the parse_vertex_type functions



class Vertex:
	def __init__(self, point, input_filename_index, projection_value, type):
		self.point = point
		self.input_filename_index = input_filename_index
		self.projection_value = projection_value
		self.type = type

	def __str__(self):
		val = str(self.point)
		val += "\n{}".format(self.type)
		# val += str(self.type)
		return val

    # vertex needs to be a type - 
    # write a class: inside vertex folder
    def is_of_type(self, type):
        """ Check if a vertex matches certain VertexType 

            :param type: A VertexType
        """
        
        return bool((self.type & type))


	# construct a vertex from a file
	# list of vertex as an object
	# has getters
	# getter: (1) is_of_type 

# data members
# {'input_filename_index': 0.0,
#  'point': array([ 0.00887632+6.3965792e-19j,  0.02405032+0.0000000e+00j,
#         -0.02531756+0.0000000e+00j]),
#  'projection_value': [(0.0004778214266161508+0j), (0.00394390844534199+0j)],
#  'type': 1}

# parse vertices

# unset = 0
# critical = 1
# semicritical = 2
# midpoint = 4
# isolated = 8
# new = 16
# curve_sample_point = 32
# surface_sample_point = 64
# removed = 128
# problematic = 256
# singular = 512