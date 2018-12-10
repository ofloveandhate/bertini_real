"""
Dan Hessler
University of Wisconsin, Eau Claire
Fall 2018
Porting to Glumpy for faster surface rendering
using OpenGL

Current Version:
    Surface agnostic
    Mouse (trackball) implementation
        Can rotate by dragging the mouse
        Can zoom with the scroll wheel
    Now with color
    Minimal implementation to change the colors given a function

TODO:
    Make the color function better
    Play with making the app interactive

"""
import numpy as np
from glumpy import app, gl, glm, gloo
from glumpy.transforms import Trackball, Position
import bertini_real

class GlumpyPlotter(object):
    """ creates the glumpyplotter object """
    def __init__(self, data=None):
        if data is None:
            self.decomposition = bertini_real.data.ReadMostRecent()
        else:
            self.decomposition = data


    def plot(self):
        """ method used to plot a surface """
        print("Plotting object of dimension: {}".format(self.decomposition.dimension))

        data = bertini_real.data.ReadMostRecent()
        tuples = data.surface.surface_sampler_data

        def extract_points(data):
            """Extract points from vertices"""
            points = []

            for vertice in data.vertices:
                """ we use 3 here because a triangle consists
                of three points """
                point = [None]*3

                for j in range(3):
                    point[j] = vertice['point'][j].real
                points.append(point)
            return points

        def make_colors(points):
            """
            computes colors according to a function!!!

            TODO:
                allow this to take in the function as a parameter
            """
            
            # f(x) = x^2 + y^2 + z^2
            colors = []

            for i in range(len(points)):

                r = points[i][0]
                g = points[i][1]
                b = points[i][2]

                # r = r**2 + g**2 + b**2
                # g = r**2 + g**2 + b**2
                # b = r**2 + g**2 + b**2

                # r = r**2 + g**2 + b**2
                # g = r + g**2 + b**2
                # b = r + g + b**2

                r = r**2
                g = g**2
                b = b**2

                colors.append([r, g, b, 1])

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

        window = app.Window(width=1024, height=1024,
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


        @window.event
        def on_draw(draw_triangles):
            """ draws the surface """
            window.clear()

            surface.draw(gl.GL_TRIANGLES, indeces)
            # surface.draw(gl.GL_LINES, indeces)

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
    """ simply calls the plot method """

    surface = GlumpyPlotter(data)
    surface.plot()
    
