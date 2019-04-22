"""
Dan Hessler
University of Wisconsin, Eau Claire
Fall 2018 - Spring 2019

.. module:: glumpyplotter
    :platform: Unix, Windows
    :synopsis: The glumpy plotter uses OpenGL to render a decomposed
               surface created using bertini_real.
"""

import math
import numpy as np
from glumpy import app, gl, glm, gloo
from glumpy.transforms import Trackball, Position
import bertini_real


class GlumpyPlotter():
    """ Creates the glumpyplotter object which is used to render
        3d surfaces created from bertini_real
    """

    def __init__(self, data=None):
        """ Reads data from disk if it is not given any.

            Args:
                data: The decomposition to be read.
        """
        if data is None:
            self.decomposition = bertini_real.data.ReadMostRecent()
        else:
            self.decomposition = data

    def plot_surface_samples(self, cmap=None, color_function=None):
        """ Method used to plot a surface. """

        data = self.decomposition
        vertex, fragment = _open_gl_code()
        sampler_data = data.surface.surface_sampler_data
        points = extract_points(data)

        triangles = np.zeros(len(sampler_data))
        for i in range(len(sampler_data)):
            for tri in sampler_data[i]:
                index_1 = int(tri[0])
                index_2 = int(tri[1])
                index_3 = int(tri[2])

                triple = [index_1, index_2, index_3]
                triangles = np.append(triangles, triple)

        window = app.Window(width=1024, height=1024,
                            color=(0.30, 0.30, 0.35, 1.00))

        verts = np.zeros(len(points), [("position", np.float32, 3),
                                       ("a_color", np.float32, 4)])
        verts["position"] = points
        verts["a_color"] = _make_colors(verts["position"],
                                        cmap,
                                        color_function)

        verts = verts.view(gloo.VertexBuffer)

        indices = np.array(triangles).astype(np.uint32)
        indices = indices.view(gloo.IndexBuffer)

        surface = gloo.Program(vertex, fragment)
        surface.bind(verts)
        surface['u_model'] = np.eye(4, dtype=np.float32)
        surface['u_view'] = glm.translation(0, 0, -5)
        surface['transform'] = Trackball(Position("position"), znear=0)
        window.attach(surface['transform'])

        @window.event
        def on_draw(draw_triangles):
            window.clear()
            surface.draw(gl.GL_TRIANGLES, indices)

        @window.event
        def on_init():
            gl.glEnable(gl.GL_DEPTH_TEST)
            gl.glPolygonOffset(1, 1)
            gl.glEnable(gl.GL_LINE_SMOOTH)
            gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)

        app.run()

    def plot_critical_curve(self, cmap=None, color_function=None):
        """ Method used to plot a surface's critical curve. """

        """
        This method does not work as intended. In its current state
        glumpy wants all of the edges to be ordered which simply can't be done.
        The idea that Danielle and I have come up with is to try and separate
        each edge with a NaN value. Instead of creating a numpy array of 0's,
        we would make a numpy array of NaN's. This was her solution for this
        problem when coding it in MATLAB. Currently, adding NaN's to the data
        causes an error -- "TypeError: object of type 'float' has no len()".
        This was caused when using the built in NaN in Python.

        Essentially, if we are able to plot all edges of a surface as separate
        entities, it should work as intended. There are examples in the glumpy
        docs which show how to plot separate objects, but implementing this
        with the current codebase seems messy...
        """

        data = self.decomposition
        vertex, fragment = _open_gl_code()

        points = extract_curve_points(data)
        critical_curve = data.surface.surface_sampler_data
        critical_curve = np.asarray(critical_curve)

        window = app.Window(width=1024, height=1024,
                            color=(0.30, 0.30, 0.35, 1.00))

        surface = gloo.Program(vertex, fragment)
        for edge in points:
            verts = np.zeros(len(edge), [("position", np.float32, 3),
                                         ("a_color", np.float32, 4)])
            verts["position"] = edge
            verts["a_color"] = _make_colors(verts["position"],
                                            cmap, color_function)
            verts = verts.view(gloo.VertexBuffer)

            # i suspect this is the problem, in that only the last one is bound
            surface.bind(verts)

        surface['u_model'] = np.eye(4, dtype=np.float32)
        surface['u_view'] = glm.translation(0, 0, -5)
        surface['transform'] = Trackball(Position("position"), znear=0)
        window.attach(surface['transform'])

        @window.event
        def on_draw(draw_triangles):
            window.clear()
            surface.draw(gl.GL_LINE_STRIP, critical_curve)

        @window.event
        def on_init():
            gl.glEnable(gl.GL_DEPTH_TEST)
            gl.glPolygonOffset(1, 1)
            gl.glEnable(gl.GL_LINE_SMOOTH)
            gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)

        app.run()


def plot_surface_samples(data=None, cmap='hsv', color_function=None):
    """
    Sets default values for colormap if none are specified.
    """
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap(cmap)

    surface = GlumpyPlotter(data)
    surface.plot_surface_samples(cmap, color_function)


def plot_critical_curve(data=None, cmap='hsv', color_function=None):
    """
    Sets default values for colormap if none are specified.
    """
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap(cmap)

    surface = GlumpyPlotter(data)
    surface.plot_critical_curve(cmap, color_function)


# --------------------------------------------------------------------------- #
#                               Helper Methods
# --------------------------------------------------------------------------- #


def default_color_function(x_coordinate, y_coordinate, z_coordinate):
    """ Helper method for make_colors()
        The default color function that is used to compute
        our points that will then be fed into the colormap.

        :param x: x coordinate of triangle.
        :param y: y coordinate of traingle.
        :param z: z coordinate of triangle.
        :rtype: A function to be used to feed the points into a colormap.

    """
    return math.sqrt(x_coordinate**2 + y_coordinate**2 + z_coordinate**2)


def _make_colors(points, cmap, color_function):
    """ Helper method for plot()
        Computes colors according to a function
        then applies a colormap from the matplotlib library.

        :param points: The triangles that will be rendered.
        :param cmap: Matplotlib colormap to be used.
        :param color_function: Function used to compute colors with colormap.
        :rtype: A list of length 4 containing R,G,B, and Alpha values.
    """
    colors = []
    data = []

    if color_function is None:
        for i in range(len(points)):
            function_result = default_color_function(points[i][0],
                                                     points[i][1],
                                                     points[i][2])
            data.append(function_result)
    else:
        for i in range(len(points)):
            function_result = color_function(points[i][0],
                                             points[i][1],
                                             points[i][2])
            data.append(function_result)

    data = np.asarray(data)

    # normalize the data
    temp = data - min(data)
    data = temp / max(temp)

    colors = cmap(data)

    return colors


def extract_points(data):
    """ Helper method for plot_surface_samples()
        Extract points from vertices

        :param data: The decomposition that we are rendering.
        :rtype: List of tuples of length 3.
    """

    points = []
    for vertex in data.vertices:
        point = [None] * 3

        for j in range(3):
            point[j] = vertex['point'][j].real
        points.append(point)

    return points


def extract_curve_points(data):
    """ Helper method for plot_critical_curve()
        Extract points from vertices

        :param data: The decomposition that we are rendering.
        :rtype: List of lists containing tuples of length 3.
    """
    # TODO turn all lists into numpy arrays, mainly for consistency and speed

    points = []
    all_points = []
    curve = []
    for edge in data.surface.critical_curve.sampler_data:
        if len(edge) <= 1:
            continue
        for vertex in edge:
            point = (data.vertices[vertex]['point'][0].real,
                     data.vertices[vertex]['point'][1].real,
                     data.vertices[vertex]['point'][2].real)
            points.append(point)
            all_points.append(point)
        points.append([np.nan, np.nan, np.nan])
        curve.append(points)
        # curve.append(float('nan'))
        points = []

    return curve


def _open_gl_code():
    # ----------------------------------------------------------------------- #
    #                    Vertex and fragment are OpenGL code
    #                             Used by Glumpy
    # ----------------------------------------------------------------------- #

    vertex = """
    attribute vec4 a_color;         // Vertex Color
    uniform mat4   u_model;         // Model matrix
    uniform mat4   u_view;          // View matrix
    uniform mat4   u_projection;    // Projection matrix
    attribute vec3 position;        // Vertex position
    varying vec4   v_color;         // Interpolated fragment color (out)
    void main()
    {
        v_color = a_color;
        gl_Position = <transform>;
    }
    """

    fragment = """
    varying vec4  v_color;          // Interpolated fragment color (in)
    void main()
    {
        gl_FragColor = v_color;
    }
    """

    return vertex, fragment
