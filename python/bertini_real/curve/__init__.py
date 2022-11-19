"""
    :platform: Unix, Windows
    :synopsis: The Curve is a child class of Decomposition. Curve contains methods that parse curve edge and samples data.
"""
import bertini_real.parse
from bertini_real.decomposition import Decomposition
import numpy as np

import copy




def is_edge_degenerate(e):
    """ 
    Check if critical point slices are degenerate (one of the endpoints is also the middle point)

        :param e: Critical point slices
        :rtype: Return True if e is degenerate
    """
    return (e[0] == e[1]) or (e[1] == e[2])

from collections import namedtuple
DirectedEdge = namedtuple('DirectedEdge', ['edge_index', 'direction'])





from enum import Enum
class EdgeDirection(Enum):
    forward = 1
    backward = 0




class CurvePiece(object):
    """
    A piece of a curve.  Must be contiguously connected.  Other people might call this a connected component.

    """
    def __init__(self, curve, edges, directed_edges):

        self.curve = curve
        self.edge_indices = edges

        # a computed / computable property.  this should be a memoized thing.
        self.directed_edges = directed_edges



    def inputfilename(self):
        return self.curve.inputfilename




    def __repr__(self):
        s = "CurvePiece on "
        s = s+f'curve {self.curve.inputfilename}, '
        s = s+f'indices {self.edge_indices}, '
        s = s+f'directed_edges {self.directed_edges}\n'
        return s


    def __str__(self):
        return repr(self)






    def to_points(self):

        print("\n\n--------\nnew to_points")
        """
        computes a list of points
        """


        if not self.directed_edges:
            raise NotImplementedError('insert code memoizing / computing the directed edges')

        # unpack a few things
        vertices = self.curve.vertices # these have already been dehomogenized
        c = self.curve

        points = np.empty((0,3))
        prev_point_index = -1

        for edge_index, direction in self.directed_edges:

            if is_edge_degenerate(c.edges[edge_index]):
                continue


            if len(c.sampler_data)>0:
                point_indices = c.sampler_data[edge_index]
            else:
                point_indices = c.edges[edge_index]

            if direction==EdgeDirection.backward:
                point_indices = point_indices[::-1] # i question the necessity of this deep copy

            print(point_indices)

            list_of_points = []
            for ii in point_indices:
                if ii != prev_point_index:
                    list_of_points.append(vertices[ii].point)
                    prev_point_index = ii


            points_this_edge = np.array(list_of_points).real
            print(points_this_edge.shape)
            print(points.shape)
            points = np.vstack((points, points_this_edge))
            print(points)
            print(f'{c.inputfilename} after vstacking',points.shape)


        return points


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


    def break_into_pieces(self, edge_indices):
        """
        

        this function returns a list of CurvePieces
        """


        print(f'\n\n-----------------\n\nbreak_into_pieces on curve "{self.inputfilename}"\n\n')
        
        import copy
        unsorted_edges = copy.deepcopy(edge_indices) # this should be a list of integers

        list_of_pieces = [] # the thing we're computing

        print(unsorted_edges)

        while len(unsorted_edges): 
            # print('\n\n================\n\ngathering data for new CurvePiece')

            directed_edges = [DirectedEdge(unsorted_edges.pop(), EdgeDirection.forward)]



            # print(edge_indices)
            # print(unsorted_edges)
            # print(directed_edges)

            added_edge = True

            while added_edge:
                # print('so far have sorted: ',directed_edges)
                added_edge = False # reset flag to keep going

                first_edge = self.edges[directed_edges[0].edge_index]
                first_edge_direction = directed_edges[0].direction
                first_point_index = first_edge[0 if first_edge_direction==EdgeDirection.forward else -1]


                last_edge = self.edges[directed_edges[-1].edge_index]
                last_edge_direction = directed_edges[-1].direction
                last_point_index = last_edge[-1 if last_edge_direction==EdgeDirection.forward else 0]

                # print(first_point_index, last_point_index, first_edge, last_edge, first_edge_direction, last_edge_direction)

                
                for edge_ind in unsorted_edges:
                    # print("testing",self.edges[edge_ind])

                    if self.edges[edge_ind][0]==last_point_index:
                        unsorted_edges.remove(edge_ind)
                        directed_edges.append(DirectedEdge(edge_ind,EdgeDirection.forward))
                        # print(f'A. added {edge_ind} going {EdgeDirection.forward}')
                        added_edge = True
                        break

                    elif self.edges[edge_ind][0]==first_point_index:
                        unsorted_edges.remove(edge_ind)
                        directed_edges.insert(0,DirectedEdge(edge_ind,EdgeDirection.backward))
                        # print(f'B. added {edge_ind} going {EdgeDirection.backward}')
                        added_edge = True
                        break

                    elif self.edges[edge_ind][-1]==last_point_index:
                        unsorted_edges.remove(edge_ind)
                        directed_edges.append(DirectedEdge(edge_ind,EdgeDirection.backward))
                        # print(f'C. added {edge_ind} going {EdgeDirection.backward}')
                        added_edge = True
                        break

                    elif self.edges[edge_ind][-1]==first_point_index:
                        unsorted_edges.remove(edge_ind)
                        directed_edges.insert(0,DirectedEdge(edge_ind,EdgeDirection.forward))
                        # print(f'D. added {edge_ind} going {EdgeDirection.forward}')
                        added_edge = True
                        break
            
            list_of_pieces.append(CurvePiece(self, [e.edge_index for e in directed_edges], directed_edges))



        print('\n\ndone sorting: \n\tunsorted')
        for edge_index in unsorted_edges:
            print(edge_index, self.edges[edge_index])

        print(f'\tsorted into {len(list_of_pieces)} pieces:')
        print(list_of_pieces)
        
        print('\n\n--------\n\n')

        return list_of_pieces



    



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
