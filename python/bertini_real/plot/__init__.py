# Nicolle Ho
# University of Notre Dame
# Spring 2017

# Silviana Amethyst
# Fall 2018, Spring 2022

# Dan Hessler
# University of Wisconsin, Eau Claire
# Fall 2018

# Foong Min Wong
# University of Wisconsin, Eau Claire
# Fall 2018

"""
Code is useful for plotting raw surfaces.  

This module is still useful, though we can now also plot surfaces using glumpy.

    :platform: Unix, Windows
    :synopsis: This module contains Plot object.
"""

import os
from bertini_real.surface import Surface, Curve, Piece
import bertini_real.util
import dill
import numpy as np
import matplotlib
# change backend with this line, if desired
# matplotlib.use('macosx')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.widgets as widgets

# print("using {} backend".format(matplotlib.get_backend()))


class StyleOptions(object):

    def __init__(self):
        self.line_thickness = 2  # there is no code using this yet.  write it.
        self.colormap = plt.cm.viridis


class VisibilityOptions(object):

    def __init__(self):
        self.vertices = False
        self.samples = False
        self.raw = False
        self.labels = False

        # for selective plotting
        self.which_faces = [] # refers to the indices in a surface or edges in a curve
        self.indices = [] #???

    def auto_adjust(self, decomposition):
        if isinstance(decomposition, Piece):
            self._adjust_for_piece(decomposition)
        elif decomposition.dimension==1:
            self._adjust_for_curve(decomposition)
        elif decomposition.dimension==2:
            self._adjust_for_surface(decomposition)
        else:
            raise NotImplementedError(f"cannot auto_adjust VisibilityOptions for dimension {decomposition.dimension} components")

    def _adjust_for_surface(self, surface):
        if len(surface.sampler_data)==0:
            self.raw = True
        else:
            self.samples = True

    def _adjust_for_piece(self, piece):
        self.which_faces = piece.indices
        self._adjust_for_surface(piece.surface)



    def _adjust_for_curve(self, curve):
        raise NotImplementedError("please implement adjusting visibility options for curves.  should be easy!")




class Options(object):

    def __init__(self):
        self.style = StyleOptions()
        self.visibility = VisibilityOptions()


class ReversableList(list):
    """ Create a ReversableList object for reversing order of data 

        :param list: The list to be read.

    """

    def reverse(self):
        return list(reversed(self))


class Plotter(object):

    def __init__(self, options=Options()):
        """ 
        Create a Plotter object for Python visualization suite 
        """

        # cache that shit, yo
        self.options = options



    def plot(self,decomposition):
        """ 
        Plot Curves/Surfaces/Pieces, axes and figures 
        """

        self.options.visibility.auto_adjust(decomposition)

        self._make_figure()
        self._make_axes(decomposition)
        self._main(decomposition)
        self._label_axes(decomposition)
        self._apply_title()

        if not self.defer_show:
            plt.show()


    def _make_widgets_surface(self,decomposition):
        """
        You must capture and store the output of this function for it to work correctly.
        """

        # first, define some actions
        def _check_actions(label):

            if label == 'Show Vertices':
                # works but with hardcoded axes
                self.options.visibility.vertices = (
                    not self.options.visibility.vertices)
                self.ax.clear()
                self._replot(decomposition)

            elif label == 'Show Smooth Surface':
                self.options.visibility.samples = (
                    not self.options.visibility.samples)
                self.ax.clear()
                self._replot(decomposition)

            elif label == 'Show Raw Surface':
                self.options.visibility.raw = (not self.options.visibility.raw)
                self.ax.clear()
                self._replot(decomposition)

            plt.draw()

        def _export_smooth_action(arg):
            if(decomposition.dimension==1):
                print('\x1b[0;31;40m'+'Unable to export OBJ file for Curve object'+'\x1b[0m')
            else:
                decomposition.export_obj_smooth()

        def _export_raw_action(arg):
            if(decomposition.dimension==1):
                print('\x1b[0;31;40m'+'Unable to export OBJ file for Curve object'+'\x1b[0m')
            else:
                decomposition.export_obj_raw()

        # measurements are in inches

        y_padding = 0.1

        button_x = 1.4
        button_y = 0.2

        check_y = 0.2 # size for ONE check
        check_x = 2.2

        inset_x = 0.1
        inset_y = 0.1

        # see https://matplotlib.org/stable/gallery/axes_grid1/demo_fixed_size_axes.html
        from mpl_toolkits.axes_grid1 import Divider, Size
        # sizes are inches

        

        
        # The width and height of the rectangle are ignored.
        

        num_check_panels = 0
        num_buttons = 0
        # First, create our check boxes

        num_checks = 3

        
        x = [Size.Fixed(inset_x), Size.Fixed(check_x)]
        y = [Size.Fixed(inset_y), Size.Fixed(check_y*num_checks)]
        divider = Divider(self.fig, (0, 0, 1, 1), x, y, aspect=False)
        check_ax = self.fig.add_axes(divider.get_position(), axes_locator=divider.new_locator(nx=1, ny=1))

        checks = widgets.CheckButtons(check_ax, ('Show Vertices', 'Show Smooth Surface', 'Show Raw Surface'),
                             (False, len(decomposition.sampler_data)>0, len(decomposition.sampler_data)==0))
        num_check_panels += 1
        checks.on_clicked(_check_actions)



        x = [Size.Fixed(inset_x), Size.Fixed(button_x)]
        y = [Size.Fixed(inset_y+(num_buttons+num_check_panels)*y_padding + check_y*num_checks + button_y*num_buttons), Size.Fixed(button_y)]
        divider = Divider(self.fig, (0, 0, 1, 1), x, y, aspect=False)
        button_smooth_ax = self.fig.add_axes(divider.get_position(), axes_locator=divider.new_locator(nx=1, ny=1))
        # button_smooth_ax = self.fig.add_axes([inset_x, inset_y+num_checks*check_y+(num_check_panels+num_buttons)*y_padding, button_x, button_y]) # as [xpos, ypos   width of sorrunding, height of surrounding]
        button_export_smooth = widgets.Button(button_smooth_ax, 'Export Smooth OBJ')
        num_buttons += 1
        button_export_smooth.on_clicked(_export_smooth_action)




        x = [Size.Fixed(inset_x), Size.Fixed(button_x)]
        y = [Size.Fixed(inset_y+(num_buttons+num_check_panels)*y_padding + check_y*num_checks + button_y*num_buttons), Size.Fixed(button_y)]
        divider = Divider(self.fig, (0, 0, 1, 1), x, y, aspect=False)
        button_raw_ax = self.fig.add_axes(divider.get_position(), axes_locator=divider.new_locator(nx=1, ny=1))

        # button_raw_ax = self.fig.add_axes([inset_x, inset_y+num_checks*check_y+num_buttons*button_y+(num_check_panels+num_buttons)*y_padding, button_x, button_y]) # as [xpos, ypos   width of sorrunding, height of surrounding]
        button_export_raw = widgets.Button(button_raw_ax, 'Export Raw OBJ')
        num_buttons += 1
        button_export_raw.on_clicked(_export_raw_action)



        return {'checks':checks, \
                'buttons':{'button_export_raw':button_export_raw, 'button_export_smooth':button_export_smooth}}

    def _main(self,decomposition):

        if isinstance(decomposition,Piece):
            self._plot_piece(decomposition)
        elif decomposition.dimension == 1:
            self._plot_curve(decomposition)
        elif decomposition.dimension == 2:
            self._plot_surface(decomposition)
        else:
            raise NotImplementedError

    def _make_figure(self):

        self.fig = plt.figure(figsize=(10,8))

    def _make_axes(self,decomposition):
        if decomposition.num_variables == 2:
            self.ax = self.fig.add_subplot(1, 1, 1)
        else:
            self.ax = self.fig.add_subplot(1, 1, 1, projection='3d')

        plt.sca(self.ax)

        try:
            self.ax.set_aspect(aspect='equal')
        except NotImplementedError as e:
            #print(e, " using `auto` instead :sadface:")
            self.ax.set_aspect(aspect='auto')


        self._adjust_axis_bounds(decomposition)

    def _adjust_axis_bounds(self,decomposition):
        d = decomposition
        self.ax.set_xlim((d.center[0]-d.radius, d.center[0]+d.radius))
        self.ax.set_ylim((d.center[1]-d.radius, d.center[1]+d.radius))

        if decomposition.num_variables == 3:
            self.ax.set_zlim((d.center[2]-d.radius, d.center[2]+d.radius))


    def _apply_title(self):
        plt.sca(self.ax)
        plt.title(os.getcwd().split(os.sep)[-1])

    def _replot(self, decomposition):
        self._main(decomposition)
        self._label_axes(decomposition)
        self._apply_title()
        self._adjust_axis_bounds(decomposition)

    def _label_axes(self, decomposition):
        # todo: these should be set from the decomposition, not assumed to be
        # x,y,z
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        if decomposition.dimension == 2:
            self.ax.set_zlabel("z")

    '''
	renders all vertices

	todo: make them colored based on a function
	'''

    def _plot_vertices(self, decomposition):
        """ Plot vertices """

        # refactored version
        xs, ys, zs = self.make_xyz(decomposition)

        if decomposition.num_variables == 2:
            verts = self.ax.scatter(xs, ys)
        else:
            verts = self.ax.scatter(xs, ys, zs, zdir='z', s=.1, alpha=1)

    # this works well for _plot_vertices
    # how can we make it work well for all the other methods???
    # returns 3 separate lists

    def make_xyz(self, decomposition):
        xs = []
        ys = []
        zs = []

        for v in decomposition.vertices:
            xs.append(v.point[0].real)
            ys.append(v.point[1].real)
            if decomposition.num_variables > 2:
                zs.append(v.point[2].real)

        return np.array(xs), np.array(ys), np.array(zs)

    def extract_points(self, decomposition):
        """ Helper method for _plot_surface_samples()
            Extract points from vertices

            :param data: Surface decomposition data
            :rtype: List of tuples of length 3.

        """
        points = []

        for vertex in decomposition.vertices:
            # allocate 3 buckets to q
            point = [None] * decomposition.num_variables

            for i in range(decomposition.num_variables):
                point[i] = vertex.point[i].real
            points.append(point)

        return points








    def _plot_curve(self, curve):
        """ 
        Plot curves 
        assumes self.options is set.  
        """
        self.points = self.extract_points(curve)
        if self.options.visibility.vertices:
            self._plot_vertices(curve)

        self._determine_nondegen_edges()

        if self.options.visibility.raw:
            self._plot_raw_edges()

        if self.options.visibility.samples:
            self._plot_edge_samples()

        

        widgets = self._make_widgets_curve(curve)


    def _plot_raw_edges(self, curve):
        """ Plot raw edges """

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
                v = curve.vertices[i]
                xs.append(v.point[0].real)
                ys.append(v.point[1].real)
                if curve.num_variables > 2:
                    zs.append(v.point[2].real)

            if curve.num_variables == 2:
                self.ax.plot(xs, ys, c=color)
            else:
                self.ax.plot(xs, ys, zs, zdir='z', c=color)

    def _plot_edge_samples(self, curve):
        """ Plot sampled edges """
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
            inds = curve.sampler_data[edge_index]
            for i in inds:
                v = curve.vertices[i]
                xs.append(v.point[0].real)
                ys.append(v.point[1].real)
                if curve.num_variables > 2:
                    zs.append(v.point[2].real)

            if curve.num_variables == 2:
                self.ax.plot(xs, ys, c=color)  # v['point'][
            else:
                self.ax.plot(xs, ys, zs, zdir='z', c=color)  # v['point']

    def _determine_nondegen_edges(self, decomposition):
        """ Determine nondegenerate edges """
        curve = decomposition
        self.nondegen = []
        for i in range(curve.num_edges):
            e = curve.edges[i]
            if e[0] != e[1] != e[2]:
                self.nondegen.append(i)




    def _plot_piece(self,piece):
        """ 
        plots a Piece of a surface.  
        """
        self.options.which_faces = piece.indices
        surf = piece.surface

        self._plot_surface(surf)



    def _plot_surface(self, surf):
        """ Plot surface"""

        self.points = self.extract_points(surf)
        if self.options.visibility.vertices:
            self._plot_vertices(surf)


        if self.options.visibility.samples:
            self._plot_surface_samples(surf)

        if self.options.visibility.raw:
            self._plot_surface_raw(surf)

        

        widgets = self._make_widgets_surface(surf)

    def _plot_surface_samples(self, surf):
        """ Plot sampler surface """
        points = self.points

        # faces = tuples
        faces = surf.sampler_data

        colormap = self.options.style.colormap

        color_list = [colormap(i) for i in np.linspace(0, 1, len(faces))]

        for i in range(len(faces)):

            color = color_list[i]

            T = []
            for tri in faces[i]:
                f = int(tri[0])
                s = int(tri[1])
                t = int(tri[2])

                k = [points[f], points[s], points[t]]

                T.append(k)

            self.ax.add_collection3d(Poly3DCollection(T, facecolors=color))
            self.ax.autoscale_view()

    def _plot_surface_raw(self, surf):
        """ Plot raw surface """
        points = self.points

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
            num_total_faces = num_total_faces + 2 * \
                (curr_face['num left'] + curr_face['num right'] + 2)
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
                # top edge
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

                    if (curr_edge[0] < 0 and curr_edge[1] < 0 and curr_edge[2] < 0):
                        continue

                    # reverse() returns None, so use ReversableList
                    curr_edge = ReversableList(curr_edge)
                    curr_edge = curr_edge.reverse()

                # bottom edge
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

                    if (curr_edge[0] < 0 and curr_edge[1] < 0 and curr_edge[2] < 0):
                        continue

                # left edge
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

                # right edge
                elif case == 4:

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

                # last case
                elif case == 5:
                    break

                # make two triangles , use the midpoint (swap the values for k)
                t1 = [points[curr_edge[0]], points[curr_edge[1]],
                      points[face['midpoint']]]
                t2 = [points[curr_edge[1]], points[curr_edge[2]],
                      points[face['midpoint']]]

                T.append(t1)
                T.append(t2)

            self.ax.add_collection3d(Poly3DCollection(T, facecolors=color))
            self.ax.autoscale_view()


def plot(data, options=Options()):
    """ Plot curve/surface

        :param data: A Curve or Surface
        :param options: style and visibility options
        :rtype: a plot
    """
    b = Plotter(options=options)

    b.plot(data)
    return b
