"""
    :platform: Unix, Windows
    :synopsis: This module contains Surface and Piece types.

"""

# Foong Min Wong
# Fall 2018 - Spring 2019
#
# Silviana Amethyst
# Spring 2022, Summer 2022
#
# University of Wisconsin, Eau Claire

import bertini_real.parse
import bertini_real.exception as br_except
import numpy as np
from bertini_real.decomposition import Decomposition
from bertini_real.curve import Curve
from bertini_real.vertex import Vertex
from bertini_real.vertex import VertexType
from bertini_real.util import ReversableList

import os
import enum
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.widgets import CheckButtons

from collections import defaultdict

import copy
import math
import trimesh




_default_file_type = 'stl'

_default_solidify_thickness =  0.1

_default_piece_basename_smooth = 'br_piece_smooth'
_default_piece_basename_raw = 'br_piece_raw'

_default_surface_basename_smooth = 'br_surface_smooth'
_default_surface_basename_raw = 'br_surface_raw'



def export_mesh(mesh, basename, autoname_using_folder=False, file_type=_default_file_type, verbose=True):
    """
    Saves a mesh (generated elsewhere) to disk,
    and returns the name of the file which was saved
    """


    if autoname_using_folder:
        fileName = os.getcwd().split(os.sep)[-1]
        outname = f'{basename}_{fileName}.{file_type}'
    else:
        outname = f'{basename}.{file_type}'

    mesh.export(file_obj=outname, file_type=file_type)

    if verbose:
        print("Exported \x1b[0;35;40m " + outname + "\x1b[0m successfully")

    return outname


def solidify_mesh(mesh, distance, offset=0):

    A = mesh # relying on referenced nature of Python here

    A.fix_normals()

    B = copy.deepcopy(mesh)

    # reverse every triangles and flip every normals
    B.invert()

    # calculate A, B vertex normals
    vertexnormsA = A.vertex_normals
    vertexnormsB = B.vertex_normals

    distA = (distance) * (offset + 1) / 2
    distB = (distance) * (1 - (offset + 1) / 2)

    # create A & B vertices that move corresponding to vertex normals and
    # distance
    A.vertices = [v + vn * distA for v,
                  vn in zip(A.vertices, A.vertex_normals)]
    B.vertices = [v + vn * distB for v,
                  vn in zip(B.vertices, B.vertex_normals)]

    numVerts = len(A.vertices)

    boundary_triangles = []

    boundary_groups = trimesh.grouping.group_rows(
        A.edges_sorted, require_count=1)

    boundary_edges = A.edges[boundary_groups]

    for edge in boundary_edges:
        for i in range(len(edge) - 1):
            t1 = [edge[i+1], edge[i], edge[i] + numVerts]
            t2 = [edge[i + 1], edge[i] + numVerts, edge[i + 1] + numVerts]

            boundary_triangles.append(t1)
            boundary_triangles.append(t2)

    Q = np.concatenate((A.vertices, B.vertices), axis=0)

    newBoundary = trimesh.Trimesh(Q, boundary_triangles)

    finalmesh = A + newBoundary + B
    return finalmesh


# i used the following to help me solve this problem:
#     https://stackoverflow.com/questions/61531935/irerate-over-package-data-files-and-copy-them-to-current-working-directory
#    lol, the misspelled "iterate" is not my fault.
#
# also:
#     https://stackoverflow.com/questions/32490629/getting-todays-date-in-yyyy-mm-dd-in-python
def copy_all_scad_files_here():
    """
    copy all source .scad files provided in bertini_real to the current directory
    """

    import pkgutil
    import pkg_resources
    from os.path import join

    scad_files = pkg_resources.resource_listdir("bertini_real", "surface/scad")

    for s in scad_files:
        contents = pkgutil.get_data('bertini_real',join('surface/scad/',s))

        with open(s,'wb') as f:
            f.write(contents)




class Piece():
    """ Create a Piece object of a surface. A surface can be made of 1 piece or multiple pieces. """

    def __init__(self, indices, surface):
        """ Initialize a Piece object with corresponding indices and surface

            :param indices: A list of nonsingular pieces' indices
            :param surface: Surface data
        """

        self.indices = indices
        self.surface = surface
        self.dimension = 2
        self.num_variables = surface.num_variables
        self.center = surface.center
        self.radius = surface.radius

    def __str__(self):
        """ toString method for Piece """
        result = "Piece of a surface with face indices:"
        result += "{}".format(self.indices)
        return result

    def __repr__(self):
        return str(self)

    def is_compact(self):
        """ Check whether a piece is:
            (1) compact (no edges touch the bounding sphere)
            (2) non-compact (at least 1 edge touches bounding sphere)


            Examples:
            sphere: (1 piece) - compact

            dingdong: (2 pieces) - one compact, one not compact

            octdong: (2 pieces) - both compact

            whitney: (2 pieces) - both non-compact
            paraboloid: (1 piece) - non compact

        """

        # bounding sphere
        sphere_curve = self.surface.sphere_curve.sampler_data

        for ii in self.indices:
            face = self.surface.faces[ii]
            if face['system top'] == "input_surf_sphere":
                return False

            if face['system bottom'] == "input_surf_sphere":
                return False

        # compact
        return True

    def centroid(self):
        """Compute the centroid of each seperate piece"""

        def flatten_and_unique(list_nD):
            """helper fucntion for centroid to get a flat list of unique values"""
            return list(set([inner for outer in list_nD for inner in outer]))

        unique_point_indices_this_piece =[]
        #face indices refer to a face on the piece
        face_indices = self.indices
        #deref the face indices to point indices and compile a list of point indices
        for i in face_indices:
             unique_point_indices_this_face=flatten_and_unique(self.surface.sampler_data[i]) #indices of points on the face
             #append the point refs on the face to the the list of all point refs on the piece
             unique_point_indices_this_piece.extend(unique_point_indices_this_face)
        #deref each point index to its point
        points = self.surface.extract_points()
        coordinates_this_piece = np.array([points[ind,:] for ind in unique_point_indices_this_piece])
        #return the mean as [x,y,z]
        return coordinates_this_piece.mean(axis=0)

    # point_singularities
    # the points on a piece ,  left and right edge will be degenerated
    # type critical

    def point_singularities(self):
        """ Compute singularity points from a Piece object

            :rtype: A list of indices of point singularities
        """

        point_singularities = []

        for face_index in self.indices:
            surf = self.surface

            # if vertex type is singular, returns true
            face = surf.faces[face_index]

            curr_edge = -10

            # top - make this to function (later)
            # top is the index, system top is the where the index lives
            if(face['system top'] == 'input_critical_curve'):
                curr_edge = surf.critical_curve.edges[face['top']]
            elif(face['system top'] == 'input_surf_sphere'):
                curr_edge = surf.sphere_curve.edges[face['top']]
            else:
                for zz in range(len(surf.singular_curves)):
                    if(surf.singular_names[zz] == face['system top']):
                        curr_edge = surf.singular_curves[zz].edges[face['top']]

            # vertices
            for ii in range(3): # 0 is left, 1 is mid, 2 is right
                if(self.surface.vertices[curr_edge[ii]].is_of_type(VertexType.singular)):
                    point_singularities.append(curr_edge[ii])


            # check bottom is singular
            if(face['system bottom'] == 'input_critical_curve'):
                curr_edge = surf.critical_curve.edges[face['bottom']]
            elif(face['system bottom'] == 'input_surf_sphere'):
                curr_edge = surf.sphere_curve.edges[face['bottom']]
            else:
                for zz in range(len(surf.singular_curves)):
                    if(surf.singular_names[zz] == face['system bottom']):
                        curr_edge = surf.singular_curves[
                            zz].edges[face['bottom']]
            # vertices
            for ii in range(3):
                if(self.surface.vertices[curr_edge[ii]].is_of_type(VertexType.singular)):
                    point_singularities.append(curr_edge[ii])

            # now we check to the left
            for edge_ind in face['left']: # this thing itself is a list
                curr_edge = surf.critical_point_slices[face['middle slice index']].edges[edge_ind]
                for ii in range(3):
                    if(self.surface.vertices[curr_edge[ii]].is_of_type(VertexType.singular)):
                        point_singularities.append(curr_edge[ii])

                        # now we check to the right
            for edge_ind in face['right']: # this thing itself is a list
                curr_edge = surf.critical_point_slices[face['middle slice index']+1].edges[edge_ind]
                for ii in range(3):
                    if(self.surface.vertices[curr_edge[ii]].is_of_type(VertexType.singular)):
                        point_singularities.append(curr_edge[ii])

        return list(set(point_singularities))

    def plot(self, color, ax):
        self.surface.plot(face_indices=self.indices, color=color, ax=ax)








    def generate_filename_no_ext(self,basename):
        """ 
        construct a filename for the piece, using face indices to make unique.  
        generates without an extension, so that it can be added later
        """

        return basename+'_'+ ('-'.join([str(i) for i in self.indices[: min(len(self.indices),3) ]])) # have to use `min` in case piece has <3 faces on it


    def generate_filename_smooth(self, file_type=_default_file_type):
        return '{}.{}'.format( self.generate_filename_no_ext(basename=_default_piece_basename_smooth), file_type)

    def generate_filename_raw(self, file_type=_default_file_type):
        return '{}.{}'.format( self.generate_filename_no_ext(basename=_default_piece_basename_raw), file_type)


    def export_smooth(self, basename=_default_piece_basename_smooth,autoname_using_folder=False,file_type=_default_file_type):

        filename_no_ext = self.generate_filename_no_ext(basename)
        self.surface.export_smooth(self.indices,filename_no_ext,autoname_using_folder,file_type)


    def export_raw(self, basename=_default_piece_basename_raw,autoname_using_folder=False,file_type=_default_file_type):

        filename_no_ext = self.generate_filename_no_ext(basename)
        self.surface.export_raw(self.indices,filename_no_ext,autoname_using_folder,file_type)




    def solidify_smooth(self, distance=_default_solidify_thickness, basename=_default_piece_basename_smooth, autoname_using_folder=False,file_type=_default_file_type):

        filename_no_ext = self.generate_filename_no_ext(basename)
        self.surface.solidify_smooth(distance, self.indices,filename_no_ext,autoname_using_folder,file_type)


    def solidify_raw(self, distance=_default_solidify_thickness, basename=_default_piece_basename_raw, autoname_using_folder=False,file_type=_default_file_type):

        filename_no_ext = self.generate_filename_no_ext(basename)
        self.surface.solidify_raw(distance, self.indices,filename_no_ext,autoname_using_folder,file_type)











class Surface(Decomposition):
    """ Create a Surface object (Child class of Decomposition)

        :param Decomposition: Decomposition data from decomp file

    """

    def __init__(self, directory, is_embedded=False,embedded_into=None):
        """ Initialize a Surface Object

            :param directory: Directory of the surface folder
        """

        self.num_faces = 0
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
        self.sampler_data = []   # store all surface_sampler data

        Decomposition.__init__(self, directory, is_embedded,embedded_into)


        # automatically parse data files to gather curve data
        self.parse_surf(self.directory)
        self.gather_faces(self.directory)
        self.gather_curves(self.directory)
        try:
            self.gather_surface_samples(self.directory)
        except:
            print("no samples found")

    def __repr__(self):
        """ toString method for Surface """
        result = "Surface with:\n"
        result += f"{self.num_faces} faces\n"
        result += f"defined using {self.num_variables} variables\n"
        result += f"center, radius of sphere: {self.center}, {self.radius}\n"
        result += f"there are {self.num_critical_slices} crit slices\n"
        result += f"and {self.num_singular_curves} singular_curves with multiplicities {self.singular_curve_multiplicities}\n\n"
        result += f"computed point set has {len(self.vertices)} total points in it\n"
        result += f""
        return result

    def __str__(self):
        return repr(self)

    def parse_surf(self, directory):
        """ Parse and store into surface data

            :param directory: Directory of the surface folder
        """
        surf_data = bertini_real.parse.parse_surf(directory)
        self.num_faces = surf_data[0]
        self.num_edges = surf_data[1]
        self.num_midpoint_slices = surf_data[2]
        self.num_critical_slices = surf_data[3]
        self.num_singular_curves = surf_data[4]
        self.singular_curve_multiplicities = surf_data[5]

    # def parse_vertex_types(self, directory):
    #     """ Parse and store vertex types data

    #     :param directory: Directory of the surface folder
    #     """
    #     vertex_types_data = bertini_real.parse.parse_vertex_types(directory)
    #     self.vertex_types_data = vertex_types_data

    def gather_faces(self, directory):
        """ Gather the faces of surface

            :param directory: Directory of the surface folder
        """
        self.faces = bertini_real.parse.parse_faces(directory)

    def gather_curves(self, directory):
        """ Gather the curves of surface

            :param directory: Directory of the surface folder
        """
        for ii in range(self.num_midpoint_slices):
            new_curve = Curve(directory + '/curve_midslice_' + str(ii),is_embedded=True,embedded_into=self)
            self.midpoint_slices.append(new_curve)
        for ii in range(self.num_critical_slices):
            new_curve = Curve(directory + '/curve_critslice_' + str(ii),is_embedded=True,embedded_into=self)
            self.critical_point_slices.append(new_curve)

        self.critical_curve = Curve(directory + '/curve_crit',is_embedded=True,embedded_into=self)
        self.sphere_curve = Curve(directory + '/curve_sphere',is_embedded=True,embedded_into=self)

        for ii in range(self.num_singular_curves):
            filename = directory + '/curve_singular_mult_' + \
                str(self.singular_curve_multiplicities[ii][0]) + '_' + str(
                    self.singular_curve_multiplicities[ii][1])
            new_curve = Curve(filename,is_embedded=True,embedded_into=self)
            self.singular_curves.append(new_curve)
            self.singular_names.append(new_curve.inputfilename)

    def gather_surface_samples(self, directory):
        """ Gather the surface samples of surface

            :param directory: Directory of the surface folder
        """
        self.sampler_data = bertini_real.parse.parse_surface_samples(
            directory)



    def check_data(self):
        """ Check data """
        try:
            if self.dimension != 2:
                print('This function designed to work on surfaces decomposed with bertini_real.  your object has dimension ' + self.dimension)

        except:
            return

    def faces_nonsingularly_connected(self, seed_index):
        """ Compute the faces that are nonsingualrly connected

            :param seed_index: Index of seed
            :rtype: Two lists containing indices of connected and unconnected faces
        """
        self.check_data()

        new_indices = [seed_index]
        connected = []

        while not(new_indices == []):
            connected.extend(new_indices)
            new_indices = self.find_connected_faces(connected)

        connected.sort()
        set_num_faces = list(range(self.num_faces))

        unconnected = list(set(set_num_faces) - set(connected))

        return connected, unconnected

    def find_connected_faces(self, current):
        """ Find connected faces from current face

            :param current: Current face
            :rtype: List containing indices of connected faces

        """

        new_indices = []

        unexamined_indices = list(range(self.num_faces))

        unexamined_indices = list(set(unexamined_indices) - set(current))

        for ii in range(len(current)):
            c = current[ii]
            f = self.faces[c]  # unpack the current face
            deleteme = []

            for jj in range(len(unexamined_indices)):
                d = unexamined_indices[jj]
                g = self.faces[d]  # unpack the examined face

                if self.faces_nonsingularly_connect(f, g):
                    new_indices.append(d)
                    deleteme.append(d)

            unexamined_indices = list(set(unexamined_indices) - set(deleteme))

        return new_indices

    def faces_nonsingularly_connect(self, f, g):
        """ Check whether faces f and g are nonsingularly connected

            :param f: Current face
            :param g: Other face
            :rtype: Return True if f and g are nonsingularly connected, else False
        """
        val = False

        if self.cannot_possibly_meet(f, g):
            return val

        elif self.faces_meet_at_left(f, g):
            val = True

        elif self.faces_meet_at_right(f, g):
            val = True

        elif self.faces_meet_at_top(f, g):
            val = True

        elif self.faces_meet_at_bottom(f, g):
            val = True

        return val

    def cannot_possibly_meet(self, f, g):
        """ Check whether faces f and g cannot possibly meet (because they are in different fiber intervals of the projection)

            :param f: Current face
            :param g: Other face
            :rtype: Return True if f and g meet, else False
        """
        val = False

        if abs(f['middle slice index'] - g['middle slice index']) >= 2:
            val = True

        return val

    def faces_meet_at_left(self, f, g):
        """ Check whether faces f and g nonsingularly connected at left

            :param f: Current face
            :param g: Other face
            :rtype: Return True if f and g nonsingularly connected at left, else False
        """
        val = False

        for ii in range(f['num left']):
            e = self.critical_point_slices[
                f['middle slice index']].edges[f['left'][ii]]
            a = e[1]

            for jj in range(g['num left']):
                E = self.critical_point_slices[
                    g['middle slice index']].edges[g['left'][jj]]
                b = E[1]

                if a == b and not(self.is_edge_degenerate(e)) and not(self.is_edge_degenerate(E)):
                    val = True
                    return val

            for jj in range(g['num right']):
                E = self.critical_point_slices[
                    g['middle slice index'] + 1].edges[g['right'][jj]]
                b = E[1]

                if a == b and not(self.is_edge_degenerate(e)) and not(self.is_edge_degenerate(E)):
                    val = True
                    return val
        return val

    def faces_meet_at_right(self, f, g):
        """ Check whether faces f and g nonsingularly connected at right

            :param f: Current face
            :param g: Other face
            :rtype: Return True if f and g nonsingularly connected at right, else False
        """
        val = False

        for ii in range(f['num right']):
            e = self.critical_point_slices[
                f['middle slice index'] + 1].edges[f['right'][ii]]
            a = e[1]

            for jj in range(g['num left']):
                E = self.critical_point_slices[
                    g['middle slice index']].edges[g['left'][jj]]
                b = E[1]

                if a == b and not(self.is_edge_degenerate(e)) and not(self.is_edge_degenerate(E)):
                    val = True
                    return val

            for jj in range(g['num right']):
                E = self.critical_point_slices[
                    g['middle slice index'] + 1].edges[g['right'][jj]]
                b = E[1]

                if a == b and not(self.is_edge_degenerate(e)) and not(self.is_edge_degenerate(E)):
                    val = True
                    return val
        return val

    def faces_meet_at_top(self, f, g):
        """ Check whether faces f and g nonsingularly connected at top

            :param f: Current face
            :param g: Other face
            :rtype: Return True if f and g nonsingularly connected at top, else False
        """
        val = False

        if(f['system top'][0:15] == 'input_singcurve'):
            return val  # cannot meet singularly, because edge is singular
        else:
            # at least they are in the same interval
            if f['middle slice index'] != g['middle slice index']:
                return val

        if (f['system top'] == g['system top']):
            if (self.critical_curve.inputfilename == f['system top']):
                if (f['top'] == g['top']):
                    val = True
                    return val

        if (f['system top'] == g['system bottom']):
            if (self.critical_curve.inputfilename == f['system top']):
                if (f['top'] == g['bottom']):
                    val = True
                    return val

        return val

    def faces_meet_at_bottom(self, f, g):
        """ Check whether faces f and g nonsingularly connected at bottom

            :param f: Current face
            :param g: Other face
            :rtype: Return True if f and g nonsingularly connected at bottom, else False
        """
        val = False

        if(f['system bottom'][0:15] == 'input_singcurve'):
            return val  # cannot meet singularly, because edge is singular
        else:
            # at least they are in the same interval
            if f['middle slice index'] != g['middle slice index']:
                return val

        if (f['system bottom'] == g['system top']):
            if (self.critical_curve.inputfilename == f['system bottom']):
                if (f['bottom'] == g['top']):
                    val = True
                    return val

        if (f['system bottom'] == g['system bottom']):
            if (self.critical_curve.inputfilename == f['system bottom']):
                if (f['bottom'] == g['bottom']):
                    val = True
                    return val

        return val

    def is_edge_degenerate(self, e):
        """ 
        Check if critical point slices are degenerate (one of the endpoints is also the middle point)

            :param e: Critical point slices
            :rtype: Return True if e is degenerate
        """
        return (e[0] == e[1]) or (e[1] == e[2])

    def separate_into_nonsingular_pieces(self):
        """ 
        Separate a surface into a list of pieces, connected at singularities
        """

        self.check_data()

        pieces = []
        connected = []
        unconnected_this = [0]

        while not(unconnected_this == []):
            seed = unconnected_this[0]
            [connected_this, unconnected_this] = self.faces_nonsingularly_connected(
                seed)
            pieces.append(Piece(connected_this, self))
            connected.extend(connected_this)
            unconnected_this = list(set(unconnected_this) - set(connected))

        return pieces

    def write_piece_data(self):
        """
        Opens and edits current scad data to set the orientation and location of a plug and socket
        """

        pieces = self.separate_into_nonsingular_pieces()

        #create a list of the centroid coordinates of each piece
        centroids = [] 
        for p in pieces:
            centroids.append(p.centroid()) 


        # compute a list of nodal singularities, and which pieces they're connected to
        pieces_connected_to_sing = defaultdict(list)
        sings_on_pieces = {} #sings are in order of the piece index

        for ii, p in enumerate(pieces):
            sing_this_piece = p.point_singularities() 
            sings_on_pieces[ii] = sing_this_piece 

            # a dictionary keyed by the integer index of the singularity, with value a list of the pieces on which it is incident
            for s in sing_this_piece: 
                pieces_connected_to_sing[s].append(ii) 


        # only put plug/socket at sing that's connected to two pieces
        #then assign that sing to k(ey) and assign [piece1,piece2] to v(value)
        wanted_sing_connections = {k:v for k,v in pieces_connected_to_sing.items() if len(v)==2}

        def unit_vector(vector):
            """Helper function to find a unit vector of a vector"""

            magnitude = np.linalg.norm(vector)
            unit=[]
            for i in range(len(vector)):
                unit.append(vector[i]/magnitude)
            return unit

        directions = defaultdict(list) # explicitly keyed by the singularities
        sing_directions = {}
        sing_locations = {}
        for sing_index,connected_pieces in wanted_sing_connections.items():
            ind_connected_piece_0 = connected_pieces[0]
            ind_connected_piece_1 = connected_pieces[1]

            # find the centroid of each piece by the index of the piece
            centroid_0 = centroids[ind_connected_piece_0]
            centroid_1 = centroids[ind_connected_piece_1]

            sing_coords =self.vertices[sing_index].point.real

            #calculate the unit vectors by traveling from the centroid of the piece to the singularity
            unit_0 = unit_vector(np.subtract(centroid_0, sing_coords))
            unit_1 = unit_vector(np.subtract(centroid_1, sing_coords))

            #find unit vector resultant of unit_0 and flipped unit_1
            direction0 = unit_vector(np.add(unit_0, np.multiply(unit_1, -1)))
            direction1 = np.multiply(direction0,-1)
            directions[sing_index] = [direction0, direction1]

            sing_directions[sing_index] = (direction0)
            sing_locations[sing_index] = (list(sing_coords))

        piece_names = []
        singularities_on_pieces = []


        singindex2int = {sing_index:ii for ii,sing_index in enumerate(wanted_sing_connections.keys())}
        int2singindex = {ii:sing_index for ii,sing_index in enumerate(wanted_sing_connections.keys())}

        # print(singindex2int)
        # print(int2singindex)

        #organize the data computed above to the scad files
        for ii,p in enumerate(pieces):
            sings_this_piece = [] 

            for sing_index,connected_pieces in wanted_sing_connections.items():
                if sing_index in sings_on_pieces[ii]: 
                    sings_this_piece.append(singindex2int[sing_index])
            singularities_on_pieces.append(sings_this_piece) 
            piece_names.append(pieces[ii].generate_filename_smooth()) 

        sing_directions_as_list = [sing_directions[sing_index] for sing_index in wanted_sing_connections.keys()]
        sing_locations_as_list = [sing_locations[sing_index] for sing_index in wanted_sing_connections.keys()]

        parity_of_sing_by_piece = [ [0 for jj in range(len(pieces))] for ii in range(len(wanted_sing_connections)) ]
        for sing_index, ps in wanted_sing_connections.items():
            parity_of_sing_by_piece[singindex2int[sing_index]][ps[0]] = -1
            parity_of_sing_by_piece[singindex2int[sing_index]][ps[1]] = 1
        
        #open and auto write the data(piece file names (without extensions), all sings of pieces, sing directions in order of sing index, sing coords in order of sing index) of the piece
        with open("br_surf_piece_data.scad", "w") as f:
            f.write(f'piece_names = [')
            f.write('"{}"'.format( '","'.join(piece_names) ))

            f.write('];\n')
            f.write(f'singularities_on_pieces = {singularities_on_pieces};\n')
            f.write(f'sing_directions = {sing_directions_as_list};\n')
            f.write(f'sing_locations = {sing_locations_as_list};\n')

            f.write(f'parities = {parity_of_sing_by_piece};\n')
            f.write(f'conn_size = 0.01;\n') #hard coded, but needs to be automatically computed

        #open up the file we just wrote to look over
        with open("br_surf_piece_data.scad", "r") as f:
            print(f.read())



    def as_mesh_smooth(self,which_faces=None):
        num_faces = self.num_faces

        if which_faces is None:
            which_faces = range(num_faces)


        points = self.extract_points()

        faces = self.sampler_data

        if not faces:
            raise br_except.SurfaceNotSampled('no surface samples found.  run sampler, or re-gather and pickle')

        vertex = []
        for p in points:
            vertex.append(p)

        vertex_np_array = np.array(vertex)

        face = []

        for ii in which_faces:
            f = faces[ii]
            for tri in f:
                face.append([tri[0], tri[1], tri[2]])

        face_np_array = np.array(face)

        A = trimesh.Trimesh(vertex_np_array, face_np_array)
        A.fix_normals()
        return A


    def as_mesh_raw(self, which_faces=None):

        num_faces = self.num_faces

        if which_faces is None:
            which_faces = range(num_faces)

        # unpack a few things
        points = self.extract_points()


        num_total_faces = 0
        for ii in range(len(which_faces)):
            curr_face = self.faces[which_faces[ii]]
            num_total_faces = num_total_faces + 2 * \
                (curr_face['num left'] + curr_face['num right'] + 2)
        num_total_faces = num_total_faces * 2

        total_face_index = 0
        TT = []


        for cc in range(len(which_faces)):
            ii = which_faces[cc]
            face = self.faces[ii]

            if (face['middle slice index']) == -1:
                continue

            case = 1
            left_edge_counter = 0
            right_edge_counter = 0

            T = []

            while 1:
                ## top edge ##
                if case == 1:

                    case += 1
                    if face['top'] < 0:
                        continue

                    curr_edge = -10
                    if(face['system top'] == 'input_critical_curve'):
                        curr_edge = self.critical_curve.edges[face['top']]
                    elif(face['system top'] == 'input_surf_sphere'):
                        curr_edge = self.sphere_curve.edges[face['top']]
                    else:
                        for zz in range(len(self.singular_curves)):
                            if(self.singular_names[zz] == face['system top']):
                                curr_edge = self.singular_curves[
                                    zz].edges[face['top']]

                    if(curr_edge == -10):
                        continue

                    if (curr_edge[0] < 0 and curr_edge[1] < 0 and curr_edge[2] < 0):
                        continue

                    curr_edge = ReversableList(curr_edge)
                    curr_edge = curr_edge.reverse()

                ## bottom edge ##
                elif case == 2:

                    case += 1

                    if face['bottom'] < 0:
                        continue

                    curr_edge = -10
                    if(face['system bottom'] == 'input_critical_curve'):
                        curr_edge = self.critical_curve.edges[face['bottom']]
                    elif(face['system bottom'] == 'input_surf_sphere'):
                        curr_edge = self.sphere_curve.edges[face['bottom']]
                    else:
                        for zz in range(len(self.singular_curves)):
                            if(self.singular_names[zz] == face['system bottom']):
                                curr_edge = self.singular_curves[
                                    zz].edges[face['bottom']]

                    if(curr_edge == -10):
                        continue

                    if (curr_edge[0] < 0 and curr_edge[1] < 0 and curr_edge[2] < 0):
                        continue

                ## left edge ##
                elif case == 3:

                    if left_edge_counter < face['num left']:

                        if face['left'][left_edge_counter] < 0:
                            continue

                        slice_ind = face['middle slice index']
                        edge_ind = face['left'][left_edge_counter]

                        curr_edge = self.critical_point_slices[
                            slice_ind].edges[edge_ind]
                        left_edge_counter = left_edge_counter + 1  # increment

                    else:
                        case = case + 1
                        continue

                ## right edge ##
                elif case == 4:

                    if right_edge_counter < face['num right']:

                        if face['right'][right_edge_counter] < 0:
                            continue

                        slice_ind = face['middle slice index'] + 1
                        edge_ind = face['right'][right_edge_counter]
                        curr_edge = self.critical_point_slices[
                            slice_ind].edges[edge_ind]
                        right_edge_counter = right_edge_counter + 1

                        curr_edge = ReversableList(curr_edge)
                        curr_edge = curr_edge.reverse()

                    else:
                        case += 1
                        continue

                ## last case ##
                elif case == 5:
                    break

                t1 = [points[curr_edge[0]], points[curr_edge[1]],
                      points[face['midpoint']]]
                t2 = [points[curr_edge[1]], points[curr_edge[2]],
                      points[face['midpoint']]]

                t3 = (curr_edge[0], curr_edge[1], face['midpoint'])
                t4 = (curr_edge[1], curr_edge[2], face['midpoint'])

                T.append(t1)
                T.append(t2)

                TT.append(t3)
                TT.append(t4)

        faces = [TT]
        vertex = []

        for p in points:
            vertex.append(p)

        vertex_np_array = np.array(vertex)
        face = []

        for f in faces:
            for tri in f:
                face.append([tri[0], tri[1], tri[2]])

        face_np_array = np.array(face)

        raw_mesh = trimesh.Trimesh(vertex_np_array, face_np_array)
        raw_mesh.fix_normals()
        return raw_mesh








    def export_raw(self, which_faces=None, basename=_default_surface_basename_raw, autoname_using_folder=False, file_type=_default_file_type):
        """ 
        Export raw decomposition of surface

        returns the name of the file which was saved
        """

        mesh = self.as_mesh_raw(which_faces)
        return export_mesh(mesh, basename, autoname_using_folder, file_type)


    def export_smooth(self, which_faces=None, basename=_default_surface_basename_smooth, autoname_using_folder=False, file_type=_default_file_type):
        """ 
        Export smooth decomposition of surface
        
        returns the name of the file which was saved
        """

        mesh = self.as_mesh_smooth(which_faces)
        return export_mesh(mesh, basename, autoname_using_folder, file_type)









    def solidify_raw(self, distance=_default_solidify_thickness, which_faces=None, basename=_default_surface_basename_raw+'solidified', autoname_using_folder=False, file_type=_default_file_type):
        """
        Solidify raw version of surface.

        Available formats include {'stl', 'obj'}.
        Default file format given by `_default_file_type`

        returns the name of the file which was saved
        """

        mesh = self.as_mesh_raw(which_faces)
        solid = solidify_mesh(mesh,distance)
        return export_mesh(solid, basename, autoname_using_folder, file_type)




    def solidify_smooth(self, distance=_default_solidify_thickness, which_faces=None, basename=_default_surface_basename_raw+'solidified', autoname_using_folder=False, file_type=_default_file_type):
        """
        Solidify smooth version of surface.  Requires that the surface has been sampled using `sampler`

        Available formats include {'stl', 'obj'}.
        Default file format given by `_default_file_type`

        returns the name of the file which was saved
        """

        mesh = self.as_mesh_smooth(which_faces)
        solid = solidify_mesh(mesh,distance)
        return export_mesh(solid, basename, autoname_using_folder, file_type)







    #################


    #  _______  _        ______     _______  _______    _______           _______  _______  _______  _______  _______
    # (  ____ \( (    /|(  __  \   (  ___  )(  ____ \  (  ____ \|\     /|(  ____ )(  ____ \(  ___  )(  ____ \(  ____ \
    # | (    \/|  \  ( || (  \  )  | (   ) || (    \/  | (    \/| )   ( || (    )|| (    \/| (   ) || (    \/| (    \/
    # | (__    |   \ | || |   ) |  | |   | || (__      | (_____ | |   | || (____)|| (__    | (___) || |      | (__
    # |  __)   | (\ \) || |   | |  | |   | ||  __)     (_____  )| |   | ||     __)|  __)   |  ___  || |      |  __)
    # | (      | | \   || |   ) |  | |   | || (              ) || |   | || (\ (   | (      | (   ) || |      | (
    # | (____/\| )  \  || (__/  )  | (___) || )        /\____) || (___) || ) \ \__| )      | )   ( || (____/\| (____/\
    # (_______/|/    )_)(______/   (_______)|/         \_______)(_______)|/   \__/|/       |/     \|(_______/(_______/



    #############
