# Nicolle Ho
# University of Notre Dame
# Spring 2017

# Danielle Brake
# Fall 2018

# Dan Hessler
# University of Wisconsin, Eau Claire
# Fall 2018

# Foong Min Wong
# University of Wisconsin, Eau Claire
# Fall 2018

import os
from bertini_real.data import BRData
from bertini_real.surface import Surface, Curve
import bertini_real.util
import dill
import numpy as np
import matplotlib
import openmesh as om
from stl import mesh
import trimesh
# change backend with this line, if desired
# matplotlib.use('macosx')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.widgets import CheckButtons
from matplotlib.widgets import Button

# print("using {} backend".format(matplotlib.get_backend()))


class StyleOptions(object):
    def __init__(self):
        self.line_thickness = 2  # there is no code using this yet.  write it.
        self.colormap = plt.cm.inferno


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
    
class Plotter(object):
    def __init__(self, data=None, options=Options()):

        if data is None:
            self.decomposition = bertini_real.data.ReadMostRecent()
        else:
            self.decomposition = data

        self.options = options

    def plot(self):

        print("plotting object of dimension "
              + str(self.decomposition.dimension))

        self.make_figure()
        self.make_axes()

        self.main()

        self.label_axes()

        # I would love to move this to its own method, but I
        # don't think that is possible with the limiations of
        # matplotlib

        # can this be done with async??

        # Create our check boxes

        # These four coordinates specify the position of the checkboxes
        rax = plt.axes([0.05, 0.4, 0.2, 0.15])
        check = CheckButtons(rax, ('Vertices', 'Surface', 'Raw Surface', 'STL'),
                             (True, True, False, False))

        # Export STL button
        # exportstl = plt.axes([0.81, 0.01, 0.15, 0.075])
        # bfvtostl = Button(exportstl,'Export STL')

        # if(bfvtostl.on_clicked()):
        #      print("Export successfully")

        def func(label):
            if label == 'Vertices':
                # works but with hardcoded axes
                self.options.visibility.vertices = (
                    not self.options.visibility.vertices)
                self.ax.clear()
                self.replot()
            elif label == 'Surface':
                self.options.visibility.samples = (
                    not self.options.visibility.samples)
                self.ax.clear()
                self.replot()
            elif label == 'Raw Surface':
                self.options.visibility.raw = (not self.options.visibility.raw)
                self.ax.clear()
                self.replot()
            elif label == 'STL':
                self.fvtostl()

            plt.draw()
        check.on_clicked(func)

        self.apply_title()

        plt.show()

    def main(self):

        self.points = self.extract_points()

        if self.options.visibility.vertices:
            self.plot_vertices()

        if self.decomposition.dimension == 1:
            self.plot_curve()
        elif self.decomposition.dimension == 2:
            self.plot_surface()

    def make_figure(self):
        self.fig = plt.figure()

    def make_axes(self):
        if self.decomposition.num_variables == 2:
            self.ax = self.fig.add_subplot(1, 1, 1)
        else:
            self.ax = self.fig.add_subplot(1, 1, 1, projection='3d')

    def apply_title(self):
        plt.title(os.getcwd().split(os.sep)[-1])

    def replot(self):
        self.main()
        self.label_axes()
        self.apply_title()

    def label_axes(self):
        # todo: these should be set from the decomposition, not assumed to be x,y,z
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        if self.decomposition.dimension is 2:
            self.ax.set_zlabel("z")

    '''
	renders all vertices

	todo: make them colored based on a function
	'''

    def plot_vertices(self):

        # refactored version
        xs, ys, zs = self.make_xyz()

        if self.decomposition.num_variables == 2:
            verts = self.ax.scatter(xs, ys)
        else:
            verts = self.ax.scatter(xs, ys, zs,
                                    zdir='z', s=.1, alpha=1)

    # this works well for plot_vertices
    # how can we make it work well for all the other methods???
    # returns 3 separate lists

    def make_xyz(self):
        xs = []
        ys = []
        zs = []

        for v in self.decomposition.vertices:
            xs.append(v['point'][0].real)
            ys.append(v['point'][1].real)
            if self.decomposition.num_variables > 2:
                zs.append(v['point'][2].real)

        return xs, ys, zs

    # would this be a better implementation than make_xyz???
    # returns one list

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

    def plot_curve(self):

        curve = self.decomposition.curve  # a local unpacking

        should_plot_raw = self.options.visibility.raw
        should_plot_samp = self.options.visibility.samples and curve.sampler_data is not None
        if self.options.visibility.autorawsamples:
            if should_plot_samp:
                should_plot_raw = False

        self.determine_nondegen_edges()

        if should_plot_raw:
            self.plot_raw_edges()

        if should_plot_samp:
            self.plot_edge_samples()

    def plot_raw_edges(self):
        curve = self.decomposition.curve  # a local unpacking

        num_nondegen_edges = len(self.nondegen)

        colormap = self.options.style.colormap
        color_list = [colormap(i)
                      for i in np.linspace(0, 1, num_nondegen_edges)]

        # instead of v['point']...etc look up into "self.points"
        for i in range(num_nondegen_edges):
            color = color_list[i]
            edge_index = self.nondegen[i]
            xs = []
            ys = []
            zs = []
            inds = curve.edges[edge_index]
            for i in inds:
                v = self.decomposition.vertices[i]
                xs.append(v['point'][0].real)
                ys.append(v['point'][1].real)
                if self.decomposition.num_variables > 2:
                    zs.append(v['point'][2].real)

            if self.decomposition.num_variables == 2:
                self.ax.plot(xs, ys, c=color)  # v['point'][
            else:
                self.ax.plot(xs, ys, zs, zdir='z', c=color)  # v['point']

    def plot_edge_samples(self):

        num_nondegen_edges = len(self.nondegen)

        colormap = self.options.style.colormap
        color_list = [colormap(i)
                      for i in np.linspace(0, 1, num_nondegen_edges)]

        for i in range(num_nondegen_edges):
            color = color_list[i]
            edge_index = self.nondegen[i]
            xs = []
            ys = []
            zs = []
            inds = self.decomposition.curve.sampler_data[edge_index]
            for i in inds:
                v = self.decomposition.vertices[i]
                xs.append(v['point'][0].real)
                ys.append(v['point'][1].real)
                if self.decomposition.num_variables > 2:
                    zs.append(v['point'][2].real)

            if self.decomposition.num_variables == 2:
                self.ax.plot(xs, ys, c=color)  # v['point'][
            else:
                self.ax.plot(xs, ys, zs, zdir='z', c=color)  # v['point']

    def determine_nondegen_edges(self):
        curve = self.decomposition.curve  # a local unpacking
        self.nondegen = []
        for i in range(curve.num_edges):
            e = curve.edges[i]
            if e[0] != e[1] != e[2]:
                self.nondegen.append(i)

    def plot_surface(self):
        surf = self.decomposition

        if self.options.visibility.samples:
            self.plot_surface_samples()

        if self.options.visibility.raw:
            self.plot_surface_raw()

    def plot_surface_samples(self):
        points = self.points

        # faces = tuples
        faces = self.decomposition.surface.surface_sampler_data

        colormap = self.options.style.colormap
        color_list = [colormap(i) for i in np.linspace(0, 1, len(faces))]
      
        for i in range(len(faces)):
            color = color_list[i]

            # Initialize T here
            T = []

            for tri in faces[i]:
                f = int(tri[0])
                s = int(tri[1])
                t = int(tri[2])

                k = [points[f], points[s], points[t]]

                T.append(k)

            self.ax.add_collection3d(Poly3DCollection(T, facecolors=color))
            
            # change this limit (resize - dynamic)
            self.ax.autoscale_view()
            # self.ax.set_xlim(-5, 5)
            # self.ax.set_ylim(-5, 5)
            # self.ax.set_zlim(-5, 5)


    def fvtostl(self):
        points = self.points
        faces = self.decomposition.surface.surface_sampler_data

        # add vertex and surface to mesh

        # mesh = om.TriMesh()

        vertex = []
        
        for p in points:
            vertex.append(p)

        vertex_np_array = np.array(vertex)

        face = []

        for f in faces: # for each face
            for tri in f: # for triangle in face
                face.append([tri[0],tri[1],tri[2]])

        face_np_array = np.array(face)
 
        obj = mesh.Mesh(np.zeros(face_np_array.shape[0], dtype=mesh.Mesh.dtype))

        for i, f in enumerate(face_np_array):
            for j in range(3):
                obj.vectors[i][j] = vertex_np_array[f[j],:]

        # fn = bertini_real.util.next_filenumber()
        fileName = os.getcwd().split(os.sep)[-1]

        obj.save('a' + fileName + '.stl')

        normmesh = trimesh.load_mesh('a' + fileName + '.stl')

        normmesh.fix_normals()

        for facet in normmesh.facets:
            normmesh.visual.face_colors[facet] = trimesh.visual.random_color()

        normmesh.show()
        
        normmesh.export(file_obj='anorm' + fileName + '.stl', file_type='stl')

        # fix_normals(normmesh,True)


        print("Export successfully")

    def plot_surface_raw(self):
        points = self.points
        surf = self.decomposition.surface
        which_faces = self.options.visibility.which_faces


        # store number of faces to num_faces
        num_faces = surf.num_faces

        # check if which_faces is empty
        if not len(which_faces):
            which_faces = list(range(num_faces))

        colormap = self.options.style.colormap
        color_list = [colormap(i) for i in np.linspace(0, 1, len(which_faces))]

        # get raw data from surface
        num_total_faces = 0
        for ii in range(len(which_faces)):
            curr_face = surf.faces[which_faces[ii]]
            num_total_faces = num_total_faces + 2 *(curr_face['num left'] + curr_face['num right'] + 2)
        num_total_faces = num_total_faces * 2
        total_face_index = 0

        for cc in range(len(which_faces)):
            color = color_list[cc]
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
                                curr_edge = surf.singular_curves[zz].edges[face['top']]
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
                                curr_edge = surf.singular_curves[zz].edges[face['bottom']]

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

                        curr_edge = surf.critical_point_slices[slice_ind].edges[edge_ind]
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
                        curr_edge = surf.critical_point_slices[slice_ind].edges[edge_ind]
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


            # while loop end heree
                # store them into objs, stl writing
                # br_plotter.fv.faces(total_face_index,:) = t1;
                # br_plotter.fv.faces(total_face_index+1,:) = t2;

                # k = [points[curr_edge[0]], points[curr_edge[1]], points[curr_edge[2]]]
                # T.append(k)
                T.append(t1)
                T.append(t2)

            self.ax.add_collection3d(Poly3DCollection(T, facecolors=color))
            self.ax.autoscale_view()
            # self.ax.set_xlim(-5, 5)
            # self.ax.set_ylim(-5, 5)
            # self.ax.set_zlim(-5, 5)


def plot(data=None, options=Options()):

    b = Plotter(data, options=options)
    b.plot()
    return b
