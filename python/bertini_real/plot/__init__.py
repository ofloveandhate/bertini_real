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
        self.linewidth = 2
        self.colormap = plt.cm.viridis
        self.colormode = ColorMode.BY_CELL

        self.mono_color = 'k'
        self.color_function = None

        self.surface_edge_color = None

        self.autotitle = True


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

        from bertini_real.vertex import VertexType
        names = [str(T).split('.')[1] for T in VertexType]

        self.vertices_by_type = {n:False for n in names}

        self.surface_samples = False
        self.surface_raw = False

        self.surface_curves = False

        self.surface_curves_raw = False
        self.surface_curves_samples = False

        surface_curve_types = ['critical','singular','midslice','critslice']
        self.surface_curves_by_type = {n:False for n in surface_curve_types}

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

        self.surface_curves = True
        self.surface_curves_by_type['singular'] = True
        self.surface_curves_by_type['critical'] = True

        if len(surface.sampler_data)==0:
            self.surface_curves_raw = True
        else:
            self.surface_curves_samples = True

        self.curve_samples = True
        if len(surface.vertices)>10000:
            print(f'have {len(surface.vertices)} vertices, so setting to invisible to start')
            self.vertices = False


    def _adjust_for_piece(self, piece):
        self._adjust_for_surface(piece.surface)



    def _adjust_for_curve(self, curve):

        if len(curve.vertices)>10000:
            print(f'have {len(curve.vertices)} vertices, so setting to invisible to start')
            self.vertices = False

        if curve.sampler_data is None:
            self.curve_raw = True
            self.curve_samples = False
        else:
            self.curve_raw = False
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

        self.surface_curves = True # turn off all curves in one stroke with this one.

        self.surface_curves_raw = True # turn off all raw curves in one stroke with this one.
        self.surface_curves_samples = True # turn off all sampled curves in one stroke with this one.

        self.surface_critical_curve = True
        self.surface_singular_curves = True
        self.surface_critical_slices = True
        self.surface_midslices = True


        # for just curves, not embedded
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

        if curve.sampler_data is None:
            self.curve_raw = True
            self.curve_samples = False
        else:
            self.curve_raw = True
            self.curve_samples = True





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


        self.init_widgets()

        self.plot_results = defaultdict(list)




    def init_widgets(self):

        self.widget_fig = None
        self.widgets = {}

        self.widgets['buttons'] = {}
        self.widgets['checks'] = {}

        self.widget_props = {}

        self.widget_props['x_padding'] = 0.1 # space between groups of items
        self.widget_props['y_padding'] = 0.1 # space between groups of items

        self.widget_props['button_x'] = 1.4
        self.widget_props['button_y'] = 0.2

        self.widget_props['check_y'] = 0.2 # size for ONE check
        self.widget_props['check_x'] = 2.2

        self.widget_props['inset_x'] = 0.1
        self.widget_props['inset_y'] = 0.1

        self.widget_props['column_x'] = 2.2 # this must be wider than any widget drawn into the columns.
        self.widget_props['column_next_y'] = defaultdict(lambda : self.widget_props['inset_y'])




    def show(self):

        self.options.render.defer_show = False
        self._adjust_all_visibility()
        plt.draw() # is this necessary???
        plt.show()

    def plot(self,decomposition):
        """ 
        Plot Curves/Surfaces/Pieces, axes and figures 
        """
        
        if self.widget_fig is None:
            self._make_new_widget_figure()


        if self.fig is None:
            self._make_new_main_figure()

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
            raise NotImplementedError("I don't know how to plot whatever you have there in that decomposition of yours.  Where did you get it?")








    #  .----..----. .---. .-. .-..----.      .--.   .---.  .---. .-. .----. .-. .-. .----.
    # { {__  | {_  {_   _}| { } || {}  }    / {} \ /  ___}{_   _}| |/  {}  \|  `| |{ {__  
    # .-._} }| {__   | |  | {_} || .--'    /  /\  \\     }  | |  | |\      /| |\  |.-._} }
    # `----' `----'  `-'  `-----'`-'       `-'  `-' `---'   `-'  `-' `----' `-' `-'`----' 

    def _make_widgets_curve(self,decomposition):


        # first, define some actions
        def _check_actions(label):

            if label == 'Vertices':
                # works but with hardcoded axes
                self.options.visibility.vertices = not self.options.visibility.vertices
                self._adjust_visibility('vertices')

            elif label == 'Smooth Curve':
                self.options.visibility.curve_samples = not self.options.visibility.curve_samples
                self._adjust_visibility('curve_samples')

            elif label == 'Raw Curve':
                self.options.visibility.curve_raw = not self.options.visibility.curve_raw
                self._adjust_visibility('curve_raw')

            if not self.options.render.defer_show:
                self.show()

        def _save_pdf(arg):
            basename = os.getcwd().split(os.sep)[-1]

            from bertini_real.util import next_filenumber
            pattern=f'{basename}*.pdf'
            n = next_filenumber(pattern)

            filename = f'{basename}{n}.pdf'
            self.fig.savefig(filename,dpi=300)

            print(f'saved with filename {filename}')



        self._add_checks_to_controller(1, 'curve_main',
            ('Smooth Curve', 'Raw Curve'),
            (decomposition.sampler_data is not None, decomposition.sampler_data is None), _check_actions)

        self._add_button_to_controller(0,'save_pdf','Save PDF',_save_pdf)



    def _make_widgets_vertices(self, decomposition):


        # otherwise, we can keep going.
        def _check_actions_vertices(vertex_type):
            # flip the bit
            self.options.visibility.vertices_by_type[vertex_type] = not self.options.visibility.vertices_by_type[vertex_type]
            
            # adjust visibility
            self._adjust_visibility_vertex_type(vertex_type)

            if not self.options.render.defer_show:
                self.show()

        def _check_actions_vertices_main(label):
            # flip the bit
            self.options.visibility.vertices = not self.options.visibility.vertices
            self._adjust_visibility('vertices')

            if not self.options.render.defer_show:
                self.show()

        names = [str(T).split('.')[1] for T in self.plot_results['vertices'].values()]
        initial_state = [self.options.visibility.vertices_by_type[n] for n in names]

        self._add_checks_to_controller(0,'vertices_by_type',names,initial_state,_check_actions_vertices)

        self._add_checks_to_controller(0, 'vertices_main',
            ('Vertices',),
            (self.options.visibility.vertices,), _check_actions_vertices_main)

    def _make_widgets_surface(self,decomposition):
        """
        You must capture and store the output of this function for it to work correctly.
        """

        # first, define some actions
        def _check_actions(label):

            if label == 'Smooth Surface':
                self.options.visibility.surface_samples = not self.options.visibility.surface_samples
                self._adjust_visibility('surface_samples')

            elif label == 'Raw Surface':
                self.options.visibility.surface_raw = not self.options.visibility.surface_raw
                self._adjust_visibility('surface_raw')

            if not self.options.render.defer_show:
                self.show()

        def _export_smooth_action(arg):
                decomposition.export_smooth()

        def _export_raw_action(arg):
                decomposition.export_raw()

        def _save_png(arg):
            basename = os.getcwd().split(os.sep)[-1]

            from bertini_real.util import next_filenumber
            pattern=f'{basename}*.png'
            n = next_filenumber(pattern)

            filename = f'{basename}{n}.png'
            self.fig.savefig(filename,dpi=300)

            print(f'saved with filename {filename}')

        

        self._add_checks_to_controller(1, 'surface main',
            ('Smooth Surface', 'Raw Surface'),
            (len(decomposition.sampler_data)>0, len(decomposition.sampler_data)==0),_check_actions)

        self._add_button_to_controller(0,'export_smooth','Export Smooth OBJ',_export_smooth_action)
        self._add_button_to_controller(0,'export_raw','Export Raw OBJ',_export_raw_action)
        self._add_button_to_controller(0,'save_png','Save PNG',_save_png)



    def _add_button_to_controller(self, column, widget_name, text, on_clicked):
        from mpl_toolkits.axes_grid1 import Divider, Size

        x = [Size.Fixed(self.widget_props['inset_x'] +column*(self.widget_props['column_x']+self.widget_props['x_padding'])), # the start position
             Size.Fixed(self.widget_props['check_x'])] # the size

        next_y = self.widget_props['column_next_y'][column]
        height = self.widget_props['button_y']

        y = [Size.Fixed(next_y), 
             Size.Fixed(height)] 

        divider = Divider(self.widget_fig, (0, 0, 1, 1), x, y, aspect=False)
        button_ax = self.widget_fig.add_axes(divider.get_position(), axes_locator=divider.new_locator(nx=1, ny=1))
        widget = widgets.Button(button_ax, text)

        widget.on_clicked(on_clicked)

        # store it so it's usable.
        self.widgets['buttons'][widget_name] = widget

        # bump so the next block in the column gets placed correctly.
        self.widget_props['column_next_y'][column] += height + self.widget_props['y_padding']
        self._resize_controller()




    def _add_checks_to_controller(self, column, widget_name, check_names, initial_values, action):
        from mpl_toolkits.axes_grid1 import Divider, Size

        assert len(check_names) == len(initial_values) # sanity check

        # compute sizes of things
        num_checks = len(check_names)

        next_y = self.widget_props['column_next_y'][column]

        x = [Size.Fixed(self.widget_props['inset_x'] +column*(self.widget_props['column_x']+self.widget_props['x_padding'])), # the start position
             Size.Fixed(self.widget_props['check_x'])] # the size

        height = self.widget_props['check_y']*num_checks + self.widget_props['y_padding']

        y = [Size.Fixed(next_y), 
             Size.Fixed(height)] 

        # make space, in a divider and axes object
        divider = Divider(self.widget_fig, (0, 0, 1, 1), x, y, aspect=False)
        ax = self.widget_fig.add_axes(divider.get_position(), axes_locator=divider.new_locator(nx=1, ny=1))

        # make the widget
        widget = widgets.CheckButtons(ax, check_names, initial_values)
        widget.on_clicked(action)

        # finally, store it
        self.widgets['checks'][widget_name] = widget

        # bump so the next block in the column gets placed correctly.
        self.widget_props['column_next_y'][column] += height + self.widget_props['y_padding']
        self._resize_controller()

    def _resize_controller(self):
        """
        makes the controller window nice and tight around the widgets we put into it :)
        """
        x = self.widget_props['column_x']*len(self.widget_props['column_next_y']) \
             + self.widget_props['x_padding']*(len(self.widget_props['column_next_y'])-1) \
             + 2*self.widget_props['inset_x']
        y = max(q for q in self.widget_props['column_next_y'].values())

        self.widget_fig.set_size_inches(x,y)



    def _make_new_widget_figure(self, figsize = (5,2)):
        """
        The default size is not very meaningful, it will be autoresized as items are put into it
        """
        self.widget_fig = plt.figure(figsize=figsize)

        self.widget_fig.canvas.mpl_connect('close_event', lambda event: plt.close(self.fig))

    def _make_new_main_figure(self, figsize = (8,8)):
        """
        One can easily override the default size by deferring showing, and setting the figure size before `plotter.show()`ing, like:

        ```
        fig = plotter.fig
        fig.set_size_inches(4, 2.75)
        ```
        """
        self.fig = plt.figure(figsize=figsize)

        self.fig.canvas.mpl_connect('close_event', lambda event: plt.close(self.widget_fig))
        print('closing one window will close both!')

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

        if self.options.style.autotitle:
            plt.title(os.getcwd().split(os.sep)[-1])

    def _adjust_all_visibility(self):
        for w in self.plot_results.keys():
            self._adjust_visibility(w)


    def _adjust_visibility(self, what):
        """
        self.show() must be called separately, otherwise get stupid results from calling this in a loop
        """

        # if what not in self.plot_results:
        #     raise RuntimeError(f"trying to adjust visibility of things in _adjust_visibility, but those things weren't rendered due to render options.  key: `{what}`.  current options: {dir(self.options.visibility)} {dir(self.options.render)}")

        if what == 'vertices':
            for T in self.plot_results['vertices'].values():
                self._adjust_visibility_vertex_type(str(T).split('.')[1])

        elif what == 'surface_curves':
            self._adjust_visibility_surface_curves(what)

        else:
            for h in self.plot_results[what]:
                h.set_visible(eval( f'self.options.visibility.{what}' ))


    def _adjust_visibility_vertex_type(self, vertex_type):
        from bertini_real.vertex import VertexType

        T = eval(f'VertexType.{vertex_type}')

        for h,t in self.plot_results['vertices'].items():  # this is a dict we're looping over
            if t == T:
                h.set_visible(self.options.visibility.vertices_by_type[vertex_type] and self.options.visibility.vertices)

    def _adjust_visibility_surface_curves(self, what):

        if not self.options.visibility.surface_curves:
            # this is the easy case, just turn them all off
            for kind,handles in self.plot_results['surface_curves'].items():
                for h in handles:
                    h.set_visible(False)

        else:
            if not self.options.visibility.surface_curves_raw:
                for kind,handles in self.plot_results['surface_curves'].items():

                    if kind.endswith('raw'):
                        for h in handles:
                            h.set_visible(False)


    def _label_axes(self, decomposition):
        # todo: these should be set from the decomposition, not assumed to be
        # x,y,z
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        if decomposition.dimension == 2:
            self.ax.set_zlabel("z")








    def _plot_vertices(self, decomposition):
        """ 
        Plot vertices 
        todo: make them colored based on a function
        """

        from bertini_real.vertex import VertexType
        import numpy as np
        
        xs, ys, zs = self.make_xyz(decomposition)

        markers = matplotlib.markers.MarkerStyle('').markers
        self.plot_results['vertices'] = {}

        for T,m in zip(VertexType, markers.keys()):

            plot_these = np.array([v.type == T for v in decomposition.vertices])

            if not np.any(plot_these):
                continue

            if decomposition.num_variables == 2:
                h = self.ax.scatter(xs[plot_these], ys[plot_these], marker=m)
            else:
                h = self.ax.scatter(xs[plot_these], ys[plot_these], zs[plot_these], zdir='z', alpha=1, marker=m)

            self.plot_results['vertices'][h] = T  # these are indexed by the handles so that they can be looped over.  i agree, it would be nice if they were flipped.

        widgets = self._make_widgets_vertices(decomposition)


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



        self._determine_nondegen_edges(curve)

        handle_name = "curve_raw"
        if self.options.render.curve_raw:
            self._plot_raw_edges(curve,handle_name)
            self._adjust_visibility(handle_name)

        handle_name = "curve_samples"
        if self.options.render.curve_samples:
            self._plot_edge_samples(curve,handle_name)
            self._adjust_visibility(handle_name)

        
        if self.options.render.vertices:
            self._plot_vertices(curve)
            self._adjust_visibility('vertices')

        widgets = self._make_widgets_curve(curve)



    def _plot_embedded_curve(self, curve, curve_name):
        """ 
        Plot an embedded curve
        assumes self.options is set.  
        """

        self.plotted_decompositions.append(curve)

        self._determine_nondegen_edges(curve)

        handle_name = curve_name+"_raw"
        if self.options.render.surface_curves_raw:
            self._plot_raw_edges(curve, curve_name=('surface_curves',handle_name))
            # self._adjust_visibility(handle_name) # for embedded, this is done at a higher level

        handle_name = curve_name+"_samples"
        if self.options.render.surface_curves_samples:
            self._plot_edge_samples(curve, curve_name=('surface_curves',handle_name))
            # self._adjust_visibility(handle_name) # for embedded, this is done at a higher level


    def _plot_raw_edges(self, curve, curve_name):
        """ Plot raw edges """


        num_nondegen_edges = len(self.nondegen)

        if self.options.style.colormode is ColorMode.BY_CELL:
            colormap = self.options.style.colormap
            color_list = [colormap(i)
                          for i in np.linspace(0, 1, num_nondegen_edges)]
        elif self.options.style.colormode is ColorMode.MONO:
            color_list = [self.options.style.mono_color]*num_nondegen_edges
        else:
            print("warning, trying to use colorfun mode for curve, but it's not implemented, thunking to mono")
            color_list = [self.options.style.mono_color]*num_nondegen_edges


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
                handle = self.ax.plot(xs, ys, c=color)
                if isinstance(curve_name,tuple):
                    self.plot_results[curve_name[0]][curve_name[1]] = handle
                else:
                    self.plot_results[curve_name].extend(handle)
            else:
                handle = self.ax.plot(xs, ys, zs, zdir='z', c=color)
                if isinstance(curve_name,tuple):
                    self.plot_results[curve_name[0]][curve_name[1]] = handle
                else:
                    self.plot_results[curve_name].extend(handle)



    def _plot_edge_samples(self, curve, curve_name):
        """ Plot sampled edges """

        num_nondegen_edges = len(self.nondegen)

        if self.options.style.colormode is ColorMode.BY_CELL:
            colormap = self.options.style.colormap
            color_list = [colormap(i)
                          for i in np.linspace(0, 1, num_nondegen_edges)]
        elif self.options.style.colormode is ColorMode.MONO:
            color_list = [self.options.style.mono_color]*num_nondegen_edges
        else:
            print("warning, trying to use colorfun mode for curve, but it's not implemented, thunking to mono")
            color_list = [self.options.style.mono_color]*num_nondegen_edges

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
                handle = self.ax.plot(xs, ys, c=color, linewidth=self.options.style.linewidth)  # v['point'][
                if isinstance(curve_name,tuple):
                    self.plot_results[curve_name[0]][curve_name[1]] = handle
                else:
                    self.plot_results[curve_name].extend(handle)
            else:
                handle = self.ax.plot(xs, ys, zs, zdir='z', c=color, linewidth=self.options.style.linewidth)  # v['point']
                if isinstance(curve_name,tuple):
                    self.plot_results[curve_name[0]][curve_name[1]] = handle
                else:
                    self.plot_results[curve_name].extend(handle)

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
        
        if not self.options.render.defer_show:
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

        if self.options.render.surface_curves:
            self._plot_surface_curves(surf)

        widgets = self._make_widgets_surface(surf)



    def _plot_surface_curves(self,surf):
        """
        plot the embedded curves in a surface
        """
        self.options.style.colormode = ColorMode.MONO

        

        if self.options.render.surface_curves:
            self.plot_results['surface_curves'] = {} # because I want to hold a dict of lists of handles

            if self.options.render.surface_critical_curve:
                self._plot_embedded_curve(surf.critical_curve, 'surface_critical_curve')

            if self.options.render.surface_singular_curves:
                for c,m in zip(surf.singular_curves, surf.singular_names):
                    self._plot_embedded_curve(c, 'surface_singular_curve')

            if self.options.render.surface_critical_slices:
                for ii,c in enumerate(surf.critical_point_slices):
                    self._plot_embedded_curve(c, f'surface_critical_slice')

            if self.options.render.surface_midslices:
                for ii,c in enumerate(surf.midpoint_slices):
                    self._plot_embedded_curve(c, f'surface_mid_slice')

            self._adjust_visibility('surface_curves')

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
            user_colors = [colormap(ii) for ii in np.linspace(0, 1, len(which_faces))]

        elif self.options.style.colormode is ColorMode.BY_FUNCTION:
            colormap = self.options.style.colormap
            color_function = self.options.style.color_function

        elif self.options.style.colormode is ColorMode.MONO:
            user_colors = [self.options.style.mono_color]*len(which_faces)

        else:
            raise NotImplementedError("unknown coloring method in style options")


        all_triangles = []
        all_colors = []

        for cc,ii in enumerate(which_faces):

            if self.options.style.colormode is ColorMode.BY_FUNCTION:
                colors_this_face = []
            else:
                colors_this_face = user_colors[cc]

            triangles_this_face = []


            for tri in faces[ii]:
                f = int(tri[0]) # i hate that these conversions are here.  this is bullshit. --sca
                s = int(tri[1])
                t = int(tri[2])

                a = points[f,:]
                b = points[s,:]
                c = points[t,:]

                triangles_this_face.append([a,b,c])

                if self.options.style.colormode is ColorMode.BY_FUNCTION:
                    colors_this_face.append( np.mean([color_function(a),color_function(b),color_function(c)], axis=0) )
            

            all_triangles.append(triangles_this_face)
            all_colors.append(colors_this_face)

        # now have all the triangles we need.  but the colors might still need some help, if using a colorfunction.
        if self.options.style.colormode is ColorMode.BY_FUNCTION:
            self._remap_colors_colorfn(all_colors)

        # finally, ready to plot
        for T,color in zip(all_triangles,all_colors):
            handle = self.ax.add_collection3d(Poly3DCollection(T, facecolors=color))
            self.plot_results['surface_samples'].append(handle)




    def _remap_colors_colorfn(self, all_colors):
        colorfunresult = all_colors[0][0]

        import numbers # https://stackoverflow.com/questions/31627321/testing-if-a-value-is-numeric

        if isinstance(colorfunresult, numbers.Number):

            # this is the 1-channel case.  it needs to get passed through the colormap.  
            # for this, the values need to be between 0 and 1 :(
            upper = max([max(c) for c in all_colors])
            lower = min([min(c) for c in all_colors])

            remap = lambda x: (x-lower)/(upper-lower)

            colormap = self.options.style.colormap

            for ii in range(len(all_colors)):
                all_colors[ii] = colormap(remap(np.array(all_colors[ii])))

        elif isinstance(colorfunresult,list) or isinstance(colorfunresult, np.ndarray):
            # this lets the user specify a function that returns 4 different values

            # first, check that we actually have 4.
            num_channels = len(colorfunresult)
            assert num_channels==4 and "your color function should return exactly 4 channels (rgba) or a scalar number (to be passed through the colormap), nothing in between"

            # do remapping to get inside 0,1
            upper = np.max([np.max(c,axis=0) for c in all_colors if c],axis=0) # the `if c` is because there might be broken faces
            lower = np.min([np.min(c,axis=0) for c in all_colors if c],axis=0)

            def deal_with_channel(l, u):
                # get the shift, scale values for a channel

                # if the min and max are the same, then we don't need to rescale, just to shift
                if u==l:
                    if l<0: return l, 1
                    if l>1: return l, 1

                    return 0, 1

                # this channel actually has span, so we need to both shift and rescale
                return l, u-l 


            subtractme = np.array([deal_with_channel(l,u)[0] for l,u in zip(lower, upper)])
            denom = np.array([deal_with_channel(l,u)[1] for l,u in zip(lower, upper)])

            # define a lambda to do the remapping into 0,1
            remap = lambda x: (x-subtractme)/denom

            # replace the colors
            for ii in range(len(all_colors)):

                if not all_colors[ii]: # to tolerate broken faces.  
                    continue

                all_colors[ii] = remap(np.array(all_colors[ii]))

        else:
            raise TypeError(f"I don't know what to do with a color function that returns an object of type {type(colorfunresult)}")


    def _plot_surface_raw(self, surf):
        """ Plot raw surface """


        # unpack a bit
        points = self.points
        which_faces = self.options.render.which_faces
        num_faces = surf.num_faces 




        # set up the colors for the faces

        if self.options.style.colormode is ColorMode.BY_CELL:
            colormap = self.options.style.colormap
            user_colors = [colormap(ii) for ii in np.linspace(0, 1, len(which_faces))]

        elif self.options.style.colormode is ColorMode.BY_FUNCTION:
            colormap = self.options.style.colormap
            color_function = self.options.style.color_function

        elif self.options.style.colormode is ColorMode.MONO:
            user_colors = [self.options.style.mono_color]*len(which_faces)

        else:
            raise NotImplementedError("unknown coloring method in style options")





        # get raw data from surface
        num_total_faces = 0
        for ii in which_faces:
            curr_face = surf.faces[ii]

            num_total_faces = num_total_faces + 2 * \
                (curr_face['num left'] + curr_face['num right'] + 2) # the last +2 is for the up/down edges, split at midpoints.
        num_total_faces = num_total_faces * 2

        all_triangles = []
        all_colors = []

        for ii in range(len(which_faces)):
            face_index = which_faces[ii]
            face = surf.faces[face_index]

            if (face['middle slice index']) == -1:
                continue
            case = 1
            left_edge_counter = 0
            right_edge_counter = 0


            triangles_this_face = []
            if self.options.style.colormode is ColorMode.BY_FUNCTION:
                colors_this_face = []
            else:
                colors_this_face = user_colors[ii]

                    
            while 1:
                # top edge
                if case == 1:
                    case += 1
                    if face['top'] < 0:
                        continue

                    curr_edge = surf.curve_with_name(face['system top']).edges[face['top']]

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

                    curr_edge = surf.curve_with_name(face['system bottom']).edges[face['bottom']]

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

                triangles_this_face.append(t1)
                triangles_this_face.append(t2)

                if self.options.style.colormode is ColorMode.BY_FUNCTION:
                    colors_this_face.append( np.mean([color_function(t1[0]),color_function(t1[1]),color_function(t1[2])], axis=0) )
                    colors_this_face.append( np.mean([color_function(t2[0]),color_function(t2[1]),color_function(t2[2])], axis=0) )

            all_triangles.append(triangles_this_face)
            all_colors.append(colors_this_face)



        if self.options.style.colormode is ColorMode.BY_FUNCTION:
            self._remap_colors_colorfn(all_colors)

        for T,color in zip(all_triangles,all_colors):
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
