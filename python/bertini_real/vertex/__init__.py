"""
.. module:: vertex
    :platform: Unix, Windows
    :synopsis: This module contains Vertex object

"""


class Vertex:
    """ Create a Vertex object 

        :param object: Curve/Surface object

    """

    def __init__(self, point, input_filename_index, projection_value, type):
        """ Initialize a Vertex object

            :param point: directory path
            :param input_filename_index
            :param projection_value: value of projection
            :param type: vertextype
        """
        self.point = point
        self.input_filename_index = input_filename_index
        self.projection_value = projection_value
        self.type = type

    def __str__(self):
        """ toString method for Vertex """
        val = str(self.point)
        val += "\n{}".format(self.type)
        # val += str(self.type)
        return val

    def is_of_type(self, type):
        """ Check if a vertex matches certain VertexType

            :param type: vertextype
        """
        # Check if a vertex matches certain VertexType
        return bool((self.type & type))
