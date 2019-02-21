"""
Foong Min Wong
University of Wisconsin, Eau Claire
Fall 2018 - Spring 2019
Implementing stereolithography (STL) export feature for Bertini_real

Current Version:


TODO:

"""
import os
import numpy as np
from bertini_real.data import BRData
from bertini_real.surface import Surface, Curve
import bertini_real
import bertini_real.util
import dill
import numpy as np
import matplotlib
from stl import mesh
import trimesh


class StyleOptions(object):

    def __init__(self):
        self.line_thickness = 2  # there is no code using this yet.  write it.


class VisibilityOptions(object):

    def __init__(self):
        self.vertices = True
        self.samples = True
        self.raw = True
        self.autorawsamples = True
        self.labels = False
        # check with matlab output
        self.which_faces = []
        self.indices = []


class Options(object):

    def __init__(self):
        self.style = StyleOptions()
        self.visibility = VisibilityOptions()


class ReversableList(list):

    def reverse(self):
        return list(reversed(self))


class NumpySTL():
    """ creates the glumpyplotter object """

    def __init__(self, data=None, options=Options()):
        # can read from the memory
        # raw data
        if data is None:
            self.decomposition = bertini_real.data.ReadMostRecent()
        else:
            self.decomposition = data

        self.options = options

    def extract_points(self):
        points = []

        # get this from plot_samples.py
        for v in self.decomposition.vertices:
            # allocate 3 buckets to q
            q = [None] * 3

            for i in range(3):
                # q[0],q[1],q[2]
                q[i] = v['point'][i].real
            points.append(q)

        return points

    def raw(self):
        print("Generating raw STL surface...")
        points = self.extract_points()

        surf = self.decomposition.surface
        which_faces = self.options.visibility.which_faces

        # store number of faces to num_faces
        num_faces = surf.num_faces

        # check if which_faces is empty
        if not len(which_faces):
            which_faces = list(range(num_faces))

        # get raw data from surface
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
                    # print('top')
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
                     # print(curr_edge)
                    if (curr_edge[0] < 0 and curr_edge[1] < 0 and curr_edge[2] < 0):
                        continue

                    # reverse() returns None, so use ReversableList
                    curr_edge = ReversableList(curr_edge)
                    curr_edge = curr_edge.reverse()

                ## bottom edge ##
                elif case == 2:
                    # print('bottom')
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

                    if (curr_edge[0] < 0 and curr_edge[1] < 0 and curr_edge[2] < 0):
                        continue

                ## left edge ##
                elif case == 3:
                    # print('left')
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
                    # print('right')
                    if right_edge_counter < face['num right']:

                        if face['right'][right_edge_counter] < 0:
                            continue

                        slice_ind = face['middle slice index'] + 1
                        edge_ind = face['right'][right_edge_counter]
                        curr_edge = surf.critical_point_slices[
                            slice_ind].edges[edge_ind]
                        right_edge_counter = right_edge_counter + 1  # increment

                        curr_edge = ReversableList(curr_edge)
                        curr_edge = curr_edge.reverse()

                    else:
                        case += 1
                        continue

                ## last case ##
                elif case == 5:
                    break

                # make two triangles , use the midpoint (swap the values for k)
                t1 = [points[curr_edge[0]], points[curr_edge[1]],
                      points[face['midpoint']]]
                t2 = [points[curr_edge[1]], points[curr_edge[2]],
                      points[face['midpoint']]]

                t3 = (curr_edge[0], curr_edge[1], face['midpoint'])
                t4 = (curr_edge[0], curr_edge[2], face['midpoint'])

                TT.append(t3)
                TT.append(t4)

                T.append(t1)
                T.append(t2)

        faces = [TT]
        vertex = []

        for p in points:
            vertex.append(p)

        vertex_np_array = np.array(vertex)
        face = []

        for f in faces:  # for each face
            for tri in f:  # for each triangle in face
                face.append([tri[0], tri[1], tri[2]])

        face_np_array = np.array(face)

        obj = mesh.Mesh(np.zeros(face_np_array.shape[
                        0], dtype=mesh.Mesh.dtype))

        for i, f in enumerate(face_np_array):
            for j in range(3):
                obj.vectors[i][j] = vertex_np_array[f[j], :]

        # get object filename
        fileName = os.getcwd().split(os.sep)[-1]

        obj.save('araw' + fileName + '.stl')

        normmesh = trimesh.load_mesh('araw' + fileName + '.stl')

        normmesh.fix_normals()

        normmesh.export(file_obj='anormraw' + fileName +
                        '.stl', file_type='stl')

        print("Export " + '\x1b[2;30;45m' + "anormraw" +
              fileName + ".stl" + '\x1b[0m' + " successfully")

    def smooth(self):
        print("Generating smooth STL surface...")
        points = self.extract_points()
        faces = self.decomposition.surface.surface_sampler_data

        # add vertex and surface to mesh
        vertex = []

        for p in points:
            vertex.append(p)

        vertex_np_array = np.array(vertex)
        face = []

        for f in faces:  # for each face
            for tri in f:  # for each triangle in face
                face.append([tri[0], tri[1], tri[2]])
                # print(tri[0], tri[1], tri[2])

        face_np_array = np.array(face)

        obj = mesh.Mesh(np.zeros(face_np_array.shape[
                        0], dtype=mesh.Mesh.dtype))

        for i, f in enumerate(face_np_array):
            for j in range(3):
                obj.vectors[i][j] = vertex_np_array[f[j], :]

        # get object filename
        fileName = os.getcwd().split(os.sep)[-1]

        obj.save('asmooth' + fileName + '.stl')

        normmesh = trimesh.load_mesh('asmooth' + fileName + '.stl')

        normmesh.fix_normals()

        for facet in normmesh.facets:
            normmesh.visual.face_colors[facet] = trimesh.visual.random_color()

        # normmesh.show()

        normmesh.export(file_obj='anormsmooth' +
                        fileName + '.stl', file_type='stl')

        print("Export " + '\x1b[2;30;45m' + "anormsmooth" +
              fileName + ".stl" + '\x1b[0m' + " successfully")


# ------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------- #

def smooth(data=None):
    surface = NumpySTL(data)
    surface.smooth()


def raw(data=None):
    surface = NumpySTL(data)
    surface.raw()
