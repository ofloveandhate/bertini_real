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


class ReversableList(list):

    def reverse(self):
        return list(reversed(self))


class GlumpyPlotter():
    """Creates the glumpyplotter object which is used to render
       3d surfaces created from bertini_real
    """

    def __init__(self, data=None):
        """Reads data from disk if it is not given any.

            Args:
                data: The decomposition to be read.
        """

        if data is None:
            self.decomposition = bertini_real.data.ReadMostRecent()
        else:
            self.decomposition = data

    def plot_surface_samples(self, cmap=None, color_function=None):
        """Method used to plot a surface."""
        data = self.decomposition

        tuples = data.surface.surface_sampler_data
        points = extract_points(data)

        triangle = []
        for i in range(len(tuples)):
            for tri in tuples[i]:
                index_1 = int(tri[0])
                index_2 = int(tri[1])
                index_3 = int(tri[2])

                triple = [index_1, index_2, index_3]
                triangle.append(triple)

        triangle = np.asarray(triangle)

        vertex, fragment = _open_gl_code()

        window = app.Window(width=2048, height=2048,
                            color=(0.30, 0.30, 0.35, 1.00))

        verts = np.zeros(len(points), [("position", np.float32, 3),
                                       ("a_color", np.float32, 4)])
        verts["position"] = points
        verts["a_color"] = make_colors(verts["position"], cmap, color_function)

        verts = verts.view(gloo.VertexBuffer)

        indices = np.array(triangle).astype(np.uint32)
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

    def plot_surface_raw(self, cmap=None, color_function=None):
        """Method used to plot a raw surface."""
        surf = self.decomposition.surface

        num_faces = surf.num_faces

        which_faces = list(range(num_faces))

        print(which_faces)

        data = self.decomposition

        tuples = data.surface.surface_sampler_data
        points = extract_points(data)

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

            while 1:
                # top edge
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

                # bottom edge
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

                # left edge
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

                # right edge
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

                # last case
                elif case == 5:
                    break

                t3 = (curr_edge[0], curr_edge[1], face['midpoint'])
                t4 = (curr_edge[1], curr_edge[2], face['midpoint'])

                TT.append(t3)
                TT.append(t4)

        triangle = []

        triangle = np.asarray(TT)

        vertex, fragment = _open_gl_code()

        window = app.Window(width=2048, height=2048,
                            color=(0.30, 0.30, 0.35, 1.00))

        verts = np.zeros(len(points), [("position", np.float32, 3),
                                       ("a_color", np.float32, 4)])
        verts["position"] = points
        verts["a_color"] = make_colors(verts["position"], cmap, color_function)

        verts = verts.view(gloo.VertexBuffer)

        indices = np.array(triangle).astype(np.uint32)
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
        """This function plots the critical curve of a surface."""
        data = self.decomposition

        tuples = data.surface.critical_curve.sampler_data
        points = extract_curve_points(data)
        triangle = tuples

        triangle = np.asarray(triangle)

        vertex, fragment = _open_gl_code()

        window = app.Window(width=2048, height=2048,
                            color=(0.30, 0.30, 0.35, 1.00))

        # DO THESE 4 LINES FOR EACH EDGE OF CRITICAL CURVE
        verts = np.zeros(len(points), [("position", np.float32, 3),
                                       ("a_color", np.float32, 4)])
        verts["position"] = points
        verts["a_color"] = make_colors(verts["position"], cmap, color_function)

        verts = verts.view(gloo.VertexBuffer)
        # stop looper here

        surface = gloo.Program(vertex, fragment)
        surface.bind(verts)
        surface['u_model'] = np.eye(4, dtype=np.float32)
        surface['u_view'] = glm.translation(0, 0, -5)
        surface['transform'] = Trackball(Position("position"), znear=0)
        window.attach(surface['transform'])

        @window.event
        def on_draw(draw_triangles):
            window.clear()
            surface.draw(gl.GL_LINE_STRIP, triangle)

        @window.event
        def on_init():
            gl.glEnable(gl.GL_DEPTH_TEST)
            gl.glPolygonOffset(1, 1)
            gl.glEnable(gl.GL_LINE_SMOOTH)
            gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)

        app.run()


# --------------------------------------------------------------------------- #

def plot_surface_samples(data=None, cmap='hsv', color_function=None):
    """Sets default values for colormap if none are specified."""
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap(cmap)

    surface = GlumpyPlotter(data)
    surface.plot_surface_samples(cmap, color_function)


def plot_surface_raw(data=None, cmap='hsv', color_function=None):
    """Sets default values for colormap if none are specified."""
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap(cmap)

    surface = GlumpyPlotter(data)
    surface.plot_surface_raw(cmap, color_function)


def plot_critical_curve(data=None, cmap='hsv', color_function=None):
    """Sets default values for colormap if none are specified."""
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


def make_colors(points, cmap, color_function):
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

    # run data through specified color function
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

    # finally, run that data through the mpl.pyplot colormap function
    # and receive our rgb values
    colors = cmap(data)

    return colors


def extract_points(data):
    """ Helper method for plot()
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
    points = []
    for edge in data.surface.critical_curve.sampler_data:
        if len(edge) <= 1:
            continue
        for vertex in edge:
            point = [None] * 3
            for i in range(3):
                point[i] = data.vertices[vertex]['point'][i].real
            points.append(point)

    return points


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
