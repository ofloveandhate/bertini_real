"""
    :platform: Unix, Windows
    :synopsis: This module contains Vertex object
"""


class Vertex:
    """ Create a Vertex object 

        :param object: Curve/Surface object

    """

    def __init__(self, point, input_filename_index, projection_value, vertex_type):
        """ Initialize a Vertex object

            :param point: coordinates
            :param input_filename_index
            :param projection_value: value of projection
            :param vertex_type: vertextype
        """
        self.point = point
        self.input_filename_index = input_filename_index
        self.projection_value = projection_value
        self.type = vertex_type

    def __repr__(self):
        """ toString method for Vertex """
        val = f"Vertex({self.point},{self.input_filename_index},{self.type})"
        return val

    def __str__(self):
        return repr(self)

    def is_of_type(self, type):
        """ Check if a vertex matches certain VertexType

            :param type: vertextype
        """
        # Check if a vertex matches certain VertexType
        return bool((self.type & type))
