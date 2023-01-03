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
from bertini_real.surface import Surface, SurfacePiece
from bertini_real.curve import Curve, CurvePiece
import bertini_real.util
from bertini_real.util import ReversableList
import dill
import numpy as np
import matplotlib
# change backend with this line, if desired
# matplotlib.use('macosx')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.widgets as widgets

# print("using {} backend".format(matplotlib.get_backend()))

from enum import Enum
from collections import defaultdict

class ColorMode(Enum):
    BY_CELL = 1
    MONO = 2
    BY_FUNCTION = 3

class StyleOptions(object):

    def __init__(self):
        self.set_defaults()

    def set_defaults(self):
        self.line_thickness = 2  # there is no code using this yet.  write it.
        self.colormap = plt.cm.viridis
        self.colormode = ColorMode.BY_CELL

        self.mono_color = None
        self.color_function = None

    def set_color_function(self, function):
        self.colormode = ColorMode.BY_FUNCTION
        self.color_function = function


class VisibilityOptions(object):
    """
    A struct-like class for serving options for the Plotter class.  
    """


    def __init__(self):
        self.set_defaults()



    def set_defaults(self):
        self.vertices = False

        self.surface_samples = False
        self.surface_raw = False

        self.curve_samples = False
        self.curve_raw = False

        self.labels = False

       

    def auto_adjust(self, decomposition):
        if isinstance(decomposition, SurfacePiece):
            self._adjust_for_piece(decomposition)
        elif isinstance(decomposition,Curve):
            self._adjust_for_curve(decomposition)
        elif isinstance(decomposition,Surface):
            self._adjust_for_surface(decomposition)
        else:
            raise NotImplementedError(f"cannot auto_adjust VisibilityOptions for dimension {decomposition.dimension} components")


        




    def _adjust_for_surface(self, surface):
        if len(surface.sampler_data)==0:
            self.surface_raw = True
        else:
            self.surface_samples = True

        if len(surface.vertices)>10000:
            self.vertices = False


    def _adjust_for_piece(self, piece):
        self._adjust_for_surface(piece.surface)



    def _adjust_for_curve(self, curve):

        if len(curve.vertices)>10000:
            self.vertices = False

        if len(curve.sampler_data)==0:
            self.curve_raw = True
        else:
            self.curve_samples = True




class RenderOptions(object):
    """
    A struct-like class for serving options for the Plotter class.  

    One curious setting is `defer_show`, which lets you stifle `show`ing results of plots.  A use of this is in case you are plotting in a loop.  
    """

    def __init__(self):
        self.set_defaults()




    def set_defaults(self):
        self.vertices = True

        self.surface_samples = True
        self.surface_raw = True

        self.curve_samples = True
        self.curve_raw = True

        self.labels = True

        # for selective plotting
        self.which_faces = [] # refers to the indices in a surface or edges in a curve
        self.which_edges = []

        self.defer_show = False


    def auto_adjust(self, decomposition):
        if isinstance(decomposition, SurfacePiece):
            self._adjust_for_piece(decomposition)
        elif isinstance(decomposition,Curve):
            self._adjust_for_curve(decomposition)
        elif isinstance(decomposition,Surface):
            self._adjust_for_surface(decomposition)
        else:
            raise NotImplementedError(f"cannot auto_adjust VisibilityOptions for dimension {decomposition.dimension} components")



    def _adjust_for_surface(self, surface):

        if len(self.which_faces)==0:
            self.which_faces = range(surface.num_faces)


    def _adjust_for_piece(self, piece):
        self.which_faces = piece.indices
        self._adjust_for_surface(piece.surface)



    def _adjust_for_curve(self, curve):
        
        if len(self.which_edges)==0:
            self.which_edges = range(curve.num_edges)





# aggregate the options.
class Options(object):

    def __init__(self):
        self.style = StyleOptions()
        self.visibility = VisibilityOptions()
        self.render = RenderOptions()







class Plotter(object):

    def __init__(self, options=Options()):
        """ 
        Create a Plotter object for Python visualization suite 
        """

        # cache that shit, yo
        self.options = options
        self.fig = None
        self.ax = None

        self.plotted_decompositions = []
        self.widgets = {}
        self.plot_results = defaultdict(list)

    def show(self):

        plt.draw() # is this necessary???
        plt.show()

    def plot(self,decomposition):
        """ 
        Plot Curves/Surfaces/Pieces, axes and figures 
        """

        if self.fig is None:
            self._make_new_figure()

        if not isinstance(decomposition,list):
            self.options.visibility.auto_adjust(decomposition)
            self.options.render.auto_adjust(decomposition)

            if self.ax is None:
                self._make_new_axes(decomposition)
                self._label_axes(decomposition)
                self._apply_title()


        self._main(decomposition)
            

        if not self.options.render.defer_show:
            self._adjust_all_visibility()
            self.show()



    def _main(self,decomposition):

        if isinstance(decomposition,list) and all([isinstance(p,SurfacePiece) for p in decomposition]):
            self._plot_pieces(decomposition)
        elif isinstance(decomposition,SurfacePiece):
            self._plot_piece(decomposition)
        elif isinstance(decomposition,Curve):
            self._plot_curve(decomposition)
        elif isinstance(decomposition,Surface):
            self._plot_surface(decomposition)
        else:
            raise NotImplementedError








    #  .----..----. .---. .-. .-..----.      .--.   .---.  .---. .-. .----. .-. .-. .----.
    # { {__  | {_  {_   _}| { } || {}  }    / {} \ /  ___}{_   _}| |/  {}  \|  `| |{ {__  
    # .-._} }| {__   | |  | {_} || .--'    /  /\  \\     }  | |  | |\      /| |\  |.-._} }
    # `----' `----'  `-'  `-----'`-'       `-'  `-' `---'   `-'  `-' `----' `-' `-'`----' 

    def _make_widgets_curve(self,decomposition):
        pass


    def _make_widgets_surface(self,decomposition):
        """
        You must capture and store the output of this function for it to work correctly.
        """

        # first, define some actions
        def _check_actions(label):

            if label == 'Vertices':
                # works but with hardcoded axes
                self.options.visibility.vertices = not self.options.visibility.vertices
                self._adjust_visibility('vertices')

            elif label == 'Smooth Surface':
                self.options.visibility.surface_samples = not self.options.visibility.surface_samples
                self._adjust_visibility('surface_samples')

            elif label == 'Raw Surface':
                self.options.visibility.surface_raw = not self.options.visibility.surface_raw
                self._adjust_visibility('surface_raw')

            self.show()

        def _export_smooth_action(arg):
            if(decomposition.dimension==1):
                raise NotImplementedError('Unable to export OBJ file for Curve object')
            else:
                decomposition.export_smooth()

        def _export_raw_action(arg):
            if(decomposition.dimension==1):
                raise NotImplementedError('Unable to export OBJ file for Curve object')
            else:
                decomposition.export_raw()

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

        
        

        num_check_panels = 0
        num_buttons = 0
        num_checks = 0
        # First, create our check boxes

        num_checks_this = 3

        x = [Size.Fixed(inset_x), Size.Fixed(check_x)]
        y = [Size.Fixed(inset_y), Size.Fixed(check_y*num_checks_this)]
        divider = Divider(self.fig, (0, 0, 1, 1), x, y, aspect=False)
        check_ax = self.fig.add_axes(divider.get_position(), axes_locator=divider.new_locator(nx=1, ny=1))

        checks = widgets.CheckButtons(check_ax, ('Vertices', 'Smooth Surface', 'Raw Surface'),
                             (False, len(decomposition.sampler_data)>0, len(decomposition.sampler_data)==0))
        num_check_panels += 1
        num_checks += num_checks_this
        checks.on_clicked(_check_actions)



        x = [Size.Fixed(inset_x), Size.Fixed(button_x)]
        y = [Size.Fixed(inset_y+(num_buttons+num_check_panels)*y_padding + check_y*num_checks + button_y*num_buttons), Size.Fixed(button_y)]
        divider = Divider(self.fig, (0, 0, 1, 1), x, y, aspect=False)
        button_smooth_ax = self.fig.add_axes(divider.get_position(), axes_locator=divider.new_locator(nx=1, ny=1))
        button_export_smooth = widgets.Button(button_smooth_ax, 'Export Smooth OBJ')
        num_buttons += 1
        button_export_smooth.on_clicked(_export_smooth_action)




        x = [Size.Fixed(inset_x), Size.Fixed(button_x)]
        y = [Size.Fixed(inset_y+(num_buttons+num_check_panels)*y_padding + check_y*num_checks + button_y*num_buttons), Size.Fixed(button_y)]
        divider = Divider(self.fig, (0, 0, 1, 1), x, y, aspect=False)
        button_raw_ax = self.fig.add_axes(divider.get_position(), axes_locator=divider.new_locator(nx=1, ny=1))
        button_export_raw = widgets.Button(button_raw_ax, 'Export Raw OBJ')
        num_buttons += 1
        button_export_raw.on_clicked(_export_raw_action)


        self.widgets['checks'] = checks
        self.widgets['buttons'] = {}
        self.widgets['buttons']['export_raw'] = button_export_raw
        self.widgets['buttons']['export_smooth'] = button_export_smooth




    def _make_new_figure(self):
        self.fig = plt.figure(figsize=(10,8))

    def _make_new_axes(self,decomposition):
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

    def _adjust_all_visibility(self):
        for w in self.plot_results.keys():
            self._adjust_visibility(w)


    def _adjust_visibility(self, what):
        """
        self.show() must be called separately, otherwise get stupid results from calling this in a loop
        """

        if what not in self.plot_results:
            raise RuntimeError(f"trying to adjust visibility of things in _adjust_visibility, but those things weren't rendered due to render options.  key: `{what}`")

        for h in self.plot_results[what]:
            h.set_visible(eval( f'self.options.visibility.{what}' ))




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
            h = self.ax.scatter(xs, ys)
        else:
            h = self.ax.scatter(xs, ys, zs, zdir='z', s=.1, alpha=1)

        self.plot_results['vertices'].append(h)

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






    ###############################
    #
    #  ,-.                 
    # /                    
    # |    . . ;-. . , ,-. 
    # \    | | |   |/  |-' 
    #  `-' `-` '   '   `-' 
                         




    def _plot_curve(self, curve):
        """ 
        Plot curves 
        assumes self.options is set.  
        """

        self.plotted_decompositions.append(curve)

        if self.options.render.vertices and not curve.is_embedded:
            self._plot_vertices(curve)

        self._determine_nondegen_edges(curve)

        if self.options.render.curve_raw:
            self._plot_raw_edges(curve)

        if self.options.render.curve_samples:
            self._plot_edge_samples(curve)

        

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




    #####################################
    #
    # ;-.                
    # |  ) o             
    # |-'  . ,-. ,-. ,-. 
    # |    | |-' |   |-' 
    # '    ' `-' `-' `-' 
                       
    def _plot_pieces(self, pieces):
        """
        A conveniece function for plotting a list of pieces.
        """

        assert( isinstance(pieces,list) and all([isinstance(p,SurfacePiece) for p in pieces]) )

        self.options.render.defer_show = True
        init_value_render_vertices = self.options.render.vertices
        self.options.render.vertices = False

        if not self.options.style.colormode is ColorMode.BY_FUNCTION:
            colormap = self.options.style.colormap
            colors = [colormap(ii) for ii in np.linspace(0, 1, len(pieces))]

            self.options.style.colormode = ColorMode.MONO



        for ii,p in enumerate(pieces):

            if not self.options.style.colormode is ColorMode.BY_FUNCTION:
                self.options.style.mono_color = colors[ii]

            self.plot(p)


        self.options.render.vertices = init_value_render_vertices

        if self.options.render.vertices:

            unique_surfaces_in_pieces = []

            for p in pieces:
                if p.surface not in unique_surfaces_in_pieces:
                    unique_surfaces_in_pieces.append(p.surface)

            for s in unique_surfaces_in_pieces: 
                self._plot_vertices(s)

        # finally, all done, so show.
        self._adjust_all_visibility()
        self.show()

    def _plot_piece(self,piece):
        """ 
        plots a Piece of a surface.  
        """

        self.plotted_decompositions.append(piece)

        self.options.render.which_faces = piece.indices
        surf = piece.surface

        self._plot_surface(surf)






    ###########################################
    #
    #  ,-.                          
    # (   `          ,-             
    #  `-.  . . ;-.  |  ,-: ,-. ,-. 
    # .   ) | | |    |- | | |   |-' 
    #  `-'  `-` '    |  `-` `-' `-' 
    #               -'              




    def _plot_surface(self, surf):
        """ 
        Plot a surface with existing options in the Plotter
        """


        self.plotted_decompositions.append(surf)

        self.points = surf.extract_points()

        if self.options.render.vertices and not surf.is_embedded:
            self._plot_vertices(surf)



        if self.options.render.surface_samples:
            self._plot_surface_samples(surf)

        if self.options.render.surface_raw:
            self._plot_surface_raw(surf)

        

        widgets = self._make_widgets_surface(surf)





    def _plot_surface_samples(self, surf):
        """ 
        Plot surface samples 
        """

        if len(surf.sampler_data)==0:
            return
            
        # locally unpack
        which_faces = self.options.render.which_faces
        points = surf.extract_points()
        faces = surf.sampler_data # these are triples of integers, indexing into the vertex_set for the decomposition



        if self.options.style.colormode is ColorMode.BY_CELL:
            colormap = self.options.style.colormap
            colors = [colormap(ii) for ii in np.linspace(0, 1, len(which_faces))]

        elif self.options.style.colormode is ColorMode.BY_FUNCTION:
            raise NotImplementedError("implement coloring by function, please")

        elif self.options.style.colormode is ColorMode.MONO:
            colors = [self.options.style.mono_color]*len(which_faces)

        else:
            raise NotImplementedError("unknown coloring method in style options")



        for cc,ii in enumerate(which_faces):

            color = colors[cc]

            T = []
            for tri in faces[ii]:
                f = int(tri[0]) # i hate that these conversions are here.  this is bullshit. --sca
                s = int(tri[1])
                t = int(tri[2])

                k = [points[f,:], points[s,:], points[t,:]]

                T.append(k)

            self.plot_results['surface_samples'].append(self.ax.add_collection3d(Poly3DCollection(T, facecolors=color)))
            # self.ax.autoscale_view()

    def _plot_surface_raw(self, surf):
        """ Plot raw surface """


        # unpack a bit
        points = self.points
        which_faces = self.options.render.which_faces
        num_faces = surf.num_faces 




        # set up the colors for the faces

        if self.options.style.colormode is ColorMode.BY_CELL:
            colormap = self.options.style.colormap
            colors = [colormap(ii) for ii in np.linspace(0, 1, len(which_faces))]

        elif self.options.style.colormode is ColorMode.BY_FUNCTION:
            raise NotImplementedError("implement coloring by function, please")

        elif self.options.style.colormode is ColorMode.MONO:
            colors = [self.options.style.mono_color]*len(which_faces)

        else:
            raise NotImplementedError("unknown coloring method in style options")





        # get raw data from surface
        num_total_faces = 0
        for ii in which_faces:
            curr_face = surf.faces[ii]

            num_total_faces = num_total_faces + 2 * \
                (curr_face['num left'] + curr_face['num right'] + 2) # the last +2 is for the up/down edges, split at midpoints.
        num_total_faces = num_total_faces * 2


        for ii in range(len(which_faces)):
            color = colors[ii]
            face_index = which_faces[ii]
            face = surf.faces[face_index]

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

            self.plot_results['surface_raw'].append(self.ax.add_collection3d(Poly3DCollection(T, facecolors=color)))
            self.ax.autoscale_view()





    #############################
    #
    # ,--.       .   ;-.  .     .   .           
    # |          |   |  ) |     |   |           
    # |-   ;-. ,-|   |-'  | ,-. |-  |-  ,-. ;-. 
    # |    | | | |   |    | | | |   |   |-' |   
    # `--' ' ' `-'   '    ' `-' `-' `-' `-' '   
    #
    ###################################



def plot(data, options=Options()):
    """ 
    Plot any of:
    * curve 
    * surface
    * piece
    * a list of pieces



        :param data: A Curve, Surface, SurfacePiece, or list of things
        :param options: style and visibility options
        :rtype: a Plotter.  
    """
    plotter = Plotter(options=options)

    plotter.plot(data)

    return plotter
