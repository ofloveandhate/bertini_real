"""
    :platform: Unix, Windows
    :synopsis: This module contains Surface and Piece types.

"""

# Foong Min Wong
# Fall 2018 - Spring 2019
#
# Silviana Amethyst
# Spring 2022
#
# University of Wisconsin, Eau Claire

import bertini_real.parse
import bertini_real.exception as br_except
import numpy as np
from bertini_real.decomposition import Decomposition
from bertini_real.curve import Curve
from bertini_real.vertex import Vertex
import bertini_real.vertextype

import os
import enum
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.widgets import CheckButtons

import copy
import math
import trimesh

class VisibilityOptions(object):

    def __init__(self):
        self.pieceVisibility = True









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
        self.visibility = VisibilityOptions()

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
        coordinates_this_piece = np.array([points[ind] for ind in unique_point_indices_this_piece])
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
                        curr_edge = surf.singular_curves[
                            zz].edges[face['top']]

            # vertices
            for ii in range(3): # 0 is left, 1 is mid, 2 is right
                if(self.surface.vertices[curr_edge[ii]].is_of_type(bertini_real.vertextype.singular)):
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
                if(self.surface.vertices[curr_edge[ii]].is_of_type(bertini_real.vertextype.singular)):
                    point_singularities.append(curr_edge[ii])

            # now we check to the left
            for edge_ind in face['left']: # this thing itself is a list
                curr_edge = surf.critical_point_slices[face['middle slice index']].edges[edge_ind]
                for ii in range(3):
                    if(self.surface.vertices[curr_edge[ii]].is_of_type(bertini_real.vertextype.singular)):
                        point_singularities.append(curr_edge[ii])

                        # now we check to the right
            for edge_ind in face['right']: # this thing itself is a list
                curr_edge = surf.critical_point_slices[face['middle slice index']+1].edges[edge_ind]
                for ii in range(3):
                    if(self.surface.vertices[curr_edge[ii]].is_of_type(bertini_real.vertextype.singular)):
                        point_singularities.append(curr_edge[ii])

        return list(set(point_singularities))

    def plot(self, color, ax):
        self.surface.plot(face_indices=self.indices, color=color, ax=ax)


    def export_obj_smooth(self, basename=f'br_piece_smooth',autoname_using_folder=False,file_type='stl'):
        if basename=='br_piece_smooth':
            basename = basename+'_'+ ('-'.join([str(i) for i in self.indices[:3]]))

        self.surface.export_obj_smooth(self.indices,basename,autoname_using_folder,file_type)

    def export_obj_raw(self, basename=f'br_piece_raw',autoname_using_folder=False,file_type='stl'):
        if basename=='br_piece_raw':
            basename = basename+'_'+ ('-'.join([str(i) for i in self.indices[:3]]))
            
        self.surface.export_obj_raw(self.indices,basename,autoname_using_folder,file_type)




############   MOVE THIS CODE
def make_figure():
    return plt.figure()

def make_axes(fig):
    return fig.add_subplot(1, 1, 1, projection='3d')

def label_axes(ax):
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

def apply_title():
    plt.title(os.getcwd().split(os.sep)[-1])

def replot(pieces, ax):
    main(pieces, ax)
    label_axes(ax)
    apply_title()

def main(pieces, ax):
    colormap = plt.cm.plasma
    # colors = plt.cm.get_cmap('hsv', len(pieces))
    color_list = [colormap(i) for i in np.linspace(0, 1, len(pieces))]

    # colors = make_colors(len(pieces))
    for ii,p in enumerate(pieces):
        if(p.visibility.pieceVisibility):
            p.plot(color=color_list[ii], ax=ax)

def plot_pieces(pieces):

    fig = make_figure()
    ax = make_axes(fig)
    label_axes(ax)

    # left, bottom, width, height
    rax = plt.axes([0.05, 0.4, 0.2, 0.05*len(pieces)])
    labels = ['piece'+str(ii) for ii,p in enumerate(pieces)]
    visibility = [True for p in enumerate(pieces)]
    check = CheckButtons(rax, labels, visibility)

    main(pieces, ax)

    def func(label):
        if label == labels[labels.index(label)]:
            pieces[labels.index(label)].visibility.pieceVisibility = (not pieces[labels.index(label)].visibility.pieceVisibility)
            ax.clear()
            replot(pieces, ax)

        plt.draw()

    check.on_clicked(func)

    apply_title()

    plt.show()

###########END MOVE THIS CODE







class Surface(Decomposition):
    """ Create a Surface object (Child class of Decomposition)

        :param Decomposition: Decomposition data from decomp file

    """

    def __init__(self, directory, is_embedded=False):
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

        Decomposition.__init__(self, directory, is_embedded)


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
            new_curve = Curve(directory + '/curve_midslice_' + str(ii),is_embedded=True)
            self.midpoint_slices.append(new_curve)
        for ii in range(self.num_critical_slices):
            new_curve = Curve(directory + '/curve_critslice_' + str(ii),is_embedded=True)
            self.critical_point_slices.append(new_curve)

        self.critical_curve = Curve(directory + '/curve_crit',is_embedded=True)
        self.sphere_curve = Curve(directory + '/curve_sphere',is_embedded=True)

        for ii in range(self.num_singular_curves):
            filename = directory + '/curve_singular_mult_' + \
                str(self.singular_curve_multiplicities[ii][0]) + '_' + str(
                    self.singular_curve_multiplicities[ii][1])
            new_curve = Curve(filename,is_embedded=True)
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

        elif self.meet_at_left(f, g):
            val = True

        elif self.meet_at_right(f, g):
            val = True

        elif self.meet_at_top(f, g):
            val = True

        elif self.meet_at_bottom(f, g):
            val = True

        return val

    def cannot_possibly_meet(self, f, g):
        """ Check whether faces f and g meet

            :param f: Current face
            :param g: Other face
            :rtype: Return True if f and g meet, else False
        """
        val = False

        if abs(f['middle slice index'] - g['middle slice index']) >= 2:
            val = True

        return val

    def meet_at_left(self, f, g):
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

                if a == b and not(self.is_degenerate(e)) and not(self.is_degenerate(E)):
                    val = True
                    return val

            for jj in range(g['num right']):
                E = self.critical_point_slices[
                    g['middle slice index'] + 1].edges[g['right'][jj]]
                b = E[1]

                if a == b and not(self.is_degenerate(e)) and not(self.is_degenerate(E)):
                    val = True
                    return val
        return val

    def meet_at_right(self, f, g):
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

                if a == b and not(self.is_degenerate(e)) and not(self.is_degenerate(E)):
                    val = True
                    return val

            for jj in range(g['num right']):
                E = self.critical_point_slices[
                    g['middle slice index'] + 1].edges[g['right'][jj]]
                b = E[1]

                if a == b and not(self.is_degenerate(e)) and not(self.is_degenerate(E)):
                    val = True
                    return val
        return val

    def meet_at_top(self, f, g):
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

    def meet_at_bottom(self, f, g):
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

    def is_degenerate(self, e):
        """ Check if critical point slices are degenerate (one of the endpoints is also the middle point)

            :param e: Critical point slices
            :rtype: Return True if e is degenerate
        """
        return (e[0] == e[1]) or (e[1] == e[2])

    def separate_into_nonsingular_pieces(self):
        """ Separate a surface into nonsingular pieces """

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

    def extract_points(self):
        """ Helper method for plot_surface_samples()
            Extract points from vertices as a list

            :param data: Surface decomposition data
            :rtype: List of tuples of length n.  This is probably terrible, and should use numpy.  Please fix this, silviana.

        """
        points = []

        for vertex in self.vertices:
            # allocate n buckets to q
            point = [None] * self.num_variables

            for i in range(self.num_variables):
                point[i] = vertex.point[i].real
            points.append(point)

        return points


    def export_obj_raw(self, which_faces=None, basename='br_surface_smooth', autoname_using_folder=False,file_type='stl'):
        return export_obj_raw(self, which_faces,basename,autoname_using_folder,file_type)

    def export_obj_smooth(self, which_faces=None, basename='br_surface_raw', autoname_using_folder=False,file_type='stl'):
        return export_obj_smooth(self, which_faces,basename,autoname_using_folder,file_type)


def separate_into_nonsingular_pieces(data=None):
    """ Separate a surface into nonsingular pieces

        :param data: Surface decomposition data. If data is None, then it reads the most recent BRData.pkl.
    """
    surface = Surface(data)
    return surface.separate_into_nonsingular_pieces()


class ReversableList(list):
    """ Create a ReversableList object for reversing order of data 

        :param list: The list to be read.

    """

    def reverse(self):
        """ Reverse function for raw surface data

            :rtype: A reversed list

        """
        return list(reversed(self))


class ObjHelper():
    """ Create a ObjHelper object for exporting OBJ files """

    def __init__(self, data=None):
        """ Read data from disk

           :param data: Surface decomposition data. If data is None, then it reads the most recent BRData.pkl.
        """

        if data is None:
            self.decomposition = bertini_real.data.read_most_recent()
        else:
            self.decomposition = data





    def export_obj_raw(self, which_faces=None, basename='br_surface_raw', autoname_using_folder=False, file_type='stl'):
        """ Export raw decomposition of surfaces to OBJ """

        print('\n' + '\x1b[0;34;40m' +
              'Generating raw surface...' + '\x1b[0m')

        points = extract_points(self)

        surf = self.decomposition

        num_faces = surf.num_faces

        if which_faces is None:
            which_faces = range(num_faces)


        num_total_faces = 0
        for ii in range(len(which_faces)):
            curr_face = surf.faces[which_faces[ii]]
            num_total_faces = num_total_faces + 2 * \
                (curr_face['num left'] + curr_face['num right'] + 2)
        num_total_faces = num_total_faces * 2

        total_face_index = 0
        TT = []


        for cc in range(len(which_faces)):
            ii = which_faces[cc]
            face = surf.faces[ii]

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
                        curr_edge = surf.critical_curve.edges[face['top']]
                    elif(face['system top'] == 'input_surf_sphere'):
                        curr_edge = surf.sphere_curve.edges[face['top']]
                    else:
                        for zz in range(len(surf.singular_curves)):
                            if(surf.singular_names[zz] == face['system top']):
                                curr_edge = surf.singular_curves[
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
                        curr_edge = surf.critical_curve.edges[face['bottom']]
                    elif(face['system bottom'] == 'input_surf_sphere'):
                        curr_edge = surf.sphere_curve.edges[face['bottom']]
                    else:
                        for zz in range(len(surf.singular_curves)):
                            if(surf.singular_names[zz] == face['system bottom']):
                                curr_edge = surf.singular_curves[
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

                        curr_edge = surf.critical_point_slices[
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
                        curr_edge = surf.critical_point_slices[
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


        if autoname_using_folder:
            fileName = os.getcwd().split(os.sep)[-1]
            outname = f'{basename}_{fileName}.{file_type}'
        else:
            outname = f'{basename}.{file_type}'


        raw_mesh.export(file_obj=outname, file_type=file_type)

        print("Exported " + '\x1b[0;35;40m' + outname + '\x1b[0m' + " successfully")




    def export_obj_smooth(self, which_faces=None, basename='br_surface_smooth', autoname_using_folder=False, file_type='stl'):
        """ Export smooth decomposition of surfaces to OBJ """

        print('\n' + '\x1b[0;34;40m' +
              'Generating smooth surface...' + '\x1b[0m')



        if which_faces is None:
            which_faces = range(num_faces)


        points = extract_points(self)

        faces = self.decomposition.sampler_data

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

        if autoname_using_folder:
            fileName = os.getcwd().split(os.sep)[-1]
            outname = f'{basename}_{fileName}.{file_type}'
        else:
            outname = f'{basename}.{file_type}'

        A.export(file_obj=outname, file_type=file_type)
        print("Exported " + '\x1b[0;35;40m' + outname + '\x1b[0m' + " successfully")





    def solidify_raw(self, which_faces=None, basename='br_surface_solidified_raw', autoname_using_folder=False, file_type='stl'):
        """ Solidify raw version of OBJ """

        print('\n' + '\x1b[0;34;40m' +
              'Solidifying raw surface...' + '\x1b[0m')

        points = extract_points(self)

        surf = self.decomposition

        num_faces = surf.num_faces


        if which_faces is None:
            which_faces = list(range(num_faces))


        num_total_faces = 0
        for ii in range(len(which_faces)):
            curr_face = surf.faces[which_faces[ii]]
            num_total_faces = num_total_faces + 2 * \
                (curr_face['num left'] + curr_face['num right'] + 2)
        num_total_faces = num_total_faces * 2

        total_face_index = 0
        TT = []

        for cc in range(len(which_faces)):
            ii = which_faces[cc]
            face = surf.faces[ii]

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
                        curr_edge = surf.critical_curve.edges[face['top']]
                    elif(face['system top'] == 'input_surf_sphere'):
                        curr_edge = surf.sphere_curve.edges[face['top']]
                    else:
                        for zz in range(len(surf.singular_curves)):
                            if(surf.singular_names[zz] == face['system top']):
                                curr_edge = surf.singular_curves[
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
                        curr_edge = surf.critical_curve.edges[face['bottom']]
                    elif(face['system bottom'] == 'input_surf_sphere'):
                        curr_edge = surf.sphere_curve.edges[face['bottom']]
                    else:
                        for zz in range(len(surf.singular_curves)):
                            if(surf.singular_names[zz] == face['system bottom']):
                                curr_edge = surf.singular_curves[
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

                        curr_edge = surf.critical_point_slices[
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
                        curr_edge = surf.critical_point_slices[
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

        A = trimesh.Trimesh(vertex_np_array, face_np_array)

        A.fix_normals()

        B = copy.deepcopy(A)

        # reverse every triangles and flip every normals
        B.invert()

        # calculate A, B vertex normals
        vertexnormsA = A.vertex_normals
        vertexnormsB = B.vertex_normals

        offset = 0
        total = 0.1

        distA = (total) * (offset + 1) / 2
        distB = (total) * (1 - (offset + 1) / 2)

        # create A & B vertices that move corresponding to vertex normals and
        # distance
        A.vertices = [v + vn * distA for v,
                      vn in zip(A.vertices, A.vertex_normals)]
        B.vertices = [v + vn * distB for v,
                      vn in zip(B.vertices, B.vertex_normals)]

        numVerts = len(A.vertices)

        T = []

        boundary_groups = trimesh.grouping.group_rows(
            A.edges_sorted, require_count=1)

        boundary_edges = A.edges[boundary_groups]

        for edge in boundary_edges:
            for i in range(len(edge) - 1):
                t1 = [edge[i], edge[i + 1], edge[i] + numVerts]
                t2 = [edge[i + 1], edge[i] + numVerts, edge[i + 1] + numVerts]

                T.append(t1)
                T.append(t2)

        Q = np.concatenate((A.vertices, B.vertices), axis=0)

        newBoundary = trimesh.Trimesh(Q, T)

        finalmesh = A + newBoundary + B

        finalmesh.fix_normals()


        if autoname_using_folder:
            fileName = os.getcwd().split(os.sep)[-1]
            outname = f'{basename}_{fileName}.{file_type}'
        else:
            outname = f'{basename}.{file_type}'

        finalmesh.export(file_obj=outname, file_type=file_type)

        print("Exported " + '\x1b[0;35;40m' + outname + '\x1b[0m' + " successfully")





    def solidify_smooth(self, which_faces=None, basename='br_surface_solidified_smooth', autoname_using_folder=False, file_type='stl'):
        """ Solidify smooth version of OBJ """

        print('\n' + '\x1b[0;34;40m' +
              'Solidifying smooth surface...' + '\x1b[0m')

        points = extract_points(self)

        faces = self.decomposition.sampler_data

        if not faces:
            raise br_except.SurfaceNotSampled('no surface samples found.  run sampler, or re-gather and pickle')



        if which_faces is None:
            which_faces = list(range(num_faces))


        vertex = []
        for p in points:
            vertex.append(p)

        vertex_np_array = np.array(vertex)

        face = []
        for f in which_faces:
            for tri in f:
                face.append([tri[0], tri[1], tri[2]])

        face_np_array = np.array(face)

        A = trimesh.Trimesh(vertex_np_array, face_np_array)

        A.fix_normals()

        B = copy.deepcopy(A)

        # reverse every triangles and flip every normals
        B.invert()

        # calculate A, B vertex normals
        vertexnormsA = A.vertex_normals
        vertexnormsB = B.vertex_normals

        offset = 0
        total = 0.1

        distA = (total) * (offset + 1) / 2
        distB = (total) * (1 - (offset + 1) / 2)

        # create A & B vertices that move corresponding to vertex normals and
        # distance
        A.vertices = [v + vn * distA for v,
                      vn in zip(A.vertices, A.vertex_normals)]
        B.vertices = [v + vn * distB for v,
                      vn in zip(B.vertices, B.vertex_normals)]

        numVerts = len(A.vertices)

        T = []

        boundary_groups = trimesh.grouping.group_rows(
            A.edges_sorted, require_count=1)

        boundary_edges = A.edges[boundary_groups]

        for edge in boundary_edges:
            for i in range(len(edge) - 1):
                t1 = [edge[i], edge[i + 1], edge[i] + numVerts]
                t2 = [edge[i + 1], edge[i] + numVerts, edge[i + 1] + numVerts]

                T.append(t1)
                T.append(t2)

        Q = np.concatenate((A.vertices, B.vertices), axis=0)

        newBoundary = trimesh.Trimesh(Q, T)

        finalmesh = A + newBoundary + B

        finalmesh.fix_normals()

        if autoname_using_folder:
            fileName = os.getcwd().split(os.sep)[-1]
            outname = f'{basename}_{fileName}.{file_type}'
        else:
            outname = f'{basename}.{file_type}'

        finalmesh.export(file_obj=outname, file_type=file_type)

        print("Exported " + '\x1b[0;35;40m' + outname + '\x1b[0m' + " successfully")





def extract_points(self):
    """ Helper method for plot_surface_samples()
        Extract points from vertices

        :param data: Surface decomposition data
        :rtype: List of tuples of length 3.
    """

    points = []

    for vertex in self.decomposition.vertices:
        point = [None] * self.decomposition.num_variables

        for i in range(self.decomposition.num_variables):
            point[i] = vertex.point[i].real
        points.append(point)

    return points


def export_obj_raw(data=None, which_faces=None, basename='br_surface_raw', autoname_using_folder=False, file_type='stl'):
    """ Export raw surface to .obj

       :param data: `Surface`, or a `Piece` of a surface. If data is `None`, then it reads the most recent `BRdataN.pkl` file from the current folder.
    """

    obj_helper = ObjHelper(data)
    obj_helper.export_obj_raw(which_faces,basename,autoname_using_folder,file_type)


def export_obj_smooth(data=None, which_faces=None, basename='br_surface_smooth', autoname_using_folder=False, file_type='stl'):
    """ Export smooth surface to .obj.  Requires that the surface has been sampled.

        :param data: `Surface`, or a `Piece`. If data is `None`, then it reads the most recent `BRdataN.pkl` file from the current folder.
    """
    obj_helper = ObjHelper(data)
    obj_helper.export_obj_smooth(which_faces,basename,autoname_using_folder,file_type)

    
def solidify_raw(data=None, which_faces=None, basename='br_surface_solidified_raw', autoname_using_folder=False, file_type='stl'):
    """ Solidify raw surface OBJ by offsetting the faces.  

        :param data: `Surface` or `Piece`. If data is `None`, then it reads the most recent `BRdataN.pkl` file from the current folder.
    """

    obj_helper = ObjHelper(data)
    obj_helper.solidify_raw(which_faces,basename,autoname_using_folder,file_type)


def solidify_smooth(data=None, which_faces=None, basename='br_surface_solidified_smooth', autoname_using_folder=False, file_type='stl'):
    """ Solidify sampled surface OBJ by offsetting the faces.  Requires that the surface has been sampled.

        :param data: `Surface` or `Piece`. If data is `None`, then it reads the most recent `BRdataN.pkl` file from the current folder.
    """

    obj_helper = ObjHelper(data)
    obj_helper.solidify_smooth(which_faces,basename,autoname_using_folder,file_type)



# def solidify(data=None, totalDist=0.1, offset=0):
#     """ Create a ObjHelper object and solidify objects """
#     surface = ObjHelper(data)
#     surface.solidify(totalDist, offset)
