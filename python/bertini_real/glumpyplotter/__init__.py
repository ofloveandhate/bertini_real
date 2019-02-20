"""
Dan Hessler
University of Wisconsin, Eau Claire
Fall 2018 - Spring 2019
Porting to Glumpy for faster surface rendering
using OpenGL

Current Version:
    Surface agnostic
    Mouse (trackball) implementation
        Can rotate by dragging the mouse
        Can zoom with the scroll wheel
    Now with color
    Now with an optional color function
        example code can be seen in the runner.py files

TODO:
    Play with making the app interactive
        i.e. checkboxes, sliders, etc
        * this would be accomplished with the pyimgui pull request on the
          glumpy repo
          This may not be possible on machines without graphics cards though
          requires OpenGL 3.2
          Macbook Pro mid 2015 has OpenGL 2.1
            * look into if this can be updated??
              not looking too promising....
"""

import numpy as np
from glumpy import app, gl, glm, gloo
from glumpy.transforms import Trackball, Position
import bertini_real
import math

class GlumpyPlotter():
    """ creates the glumpyplotter object """
    def __init__(self, data=None):
        if data is None:
            self.decomposition = bertini_real.data.ReadMostRecent()
        else:
            self.decomposition = data


    def plot(self):
        """ method used to plot a surface """
        print("Plotting object of dimension: {}".format(self.decomposition.dimension))

        data = self.decomposition
        tuples = data.surface.surface_sampler_data

        def extract_points(data):
            """Extract points from vertices"""

            points = []

            for vertex in data.vertices:
                """ we use 3 here because a triangle consists
                of three points """
                point = [None]*3

                for j in range(3):
                    point[j] = vertex['point'][j].real
                points.append(point)
            return points

        #  def make_colors(points, fn = default_color_function()):
        def make_colors(points):
            """
            computes colors according to a function
            then applies a colormap from the matplotlib library
            """
            import matplotlib.pyplot as plt

            colors = []

            # first, loop through all of the points and get a list
            # of values returned from this function
            def default_color_function(x,y,z):
                return math.sqrt(x**2 + y**2 + z**2)

            data = []
            for i in range(len(points)):
                function_result = default_color_function(points[i][0],points[i][1],points[i][2])
                data.append(function_result)
            data = np.asarray(data)

            # next, normalize the data
            """
            looks like its already normalized?? will come back to this later
            """

            # soething like
            temp1 = data-min(data)
            data = temp1 / max(temp1)
            # print(data)

            # finally, run that data through the mpl.pyplot colormap function
            # and receive our rgb values
            cmap = plt.cm.hsv
            colors = cmap(data)

            return colors

        points = extract_points(data)

        triangle = []

        for i in range(len(tuples)):
            for tri in tuples[i]:
                f = int(tri[0])
                s = int(tri[1])
                t = int(tri[2])

                k = [f, s, t]
                triangle.append(k)

        triangle = np.asarray(triangle)

# ------------------------------------------------------------------------------------- #
# Vertex and fragment are OpenGL code

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
            // gl_FragColor = vec4(1.0, 1.0, 0.0, 1.0);
            gl_FragColor = v_color;
        }
        """

        window = app.Window(width=2048, height=2048,
                            color=(0.30, 0.30, 0.35, 1.00))


        verts = np.zeros(len(points), [("position", np.float32, 3),
                                        ("a_color", np.float32, 4)])
        verts["position"] = points
        verts["a_color"] = make_colors(verts["position"])

        verts = verts.view(gloo.VertexBuffer)
        indeces = np.array(triangle).astype(np.uint32)
        indeces = indeces.view(gloo.IndexBuffer)

        surface = gloo.Program(vertex, fragment)
        surface.bind(verts)
        surface['u_model'] = np.eye(4, dtype=np.float32)
        surface['u_view'] = glm.translation(0, 0, -5)
        surface['transform'] = Trackball(Position("position"))
        window.attach(surface['transform'])
        framebuffer = gloo.FrameBuffer()


        @window.event
        def on_draw(draw_triangles):
            """ draws the surface """
            window.clear()

            #  framebuffer.activate()
            surface.draw(gl.GL_TRIANGLES, indeces)
            #  surface.draw(gl.GL_LINES, indeces)
            #  framebuffer.deactivate()

        @window.event
        def on_init():
            """ settings for OpenGL, not sure what they all do """
            gl.glEnable(gl.GL_DEPTH_TEST)
            gl.glPolygonOffset(1, 1)
            gl.glEnable(gl.GL_LINE_SMOOTH)
            gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)


        #  app.run(interactive=True)
        app.run()

# ------------------------------------------------------------------------------------- #

def plot(data=None):
    """
    simply calls the plot method
    color_function contains 3 functions to compute the colors
    of the surface
    """
    surface = GlumpyPlotter(data)
    surface.plot()
