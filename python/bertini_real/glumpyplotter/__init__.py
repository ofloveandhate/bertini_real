"""
Dan Hessler
University of Wisconsin, Eau Claire
Fall 2018
Porting to Glumpy for faster surface rendering
using OpenGL
"""
import numpy as np
from glumpy import app, gl, glm, gloo
from glumpy.transforms import Trackball, Position
import bertini_real

THETA, PHI = 40, 30


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

        vertex = """
        uniform mat4   u_model;         // Model matrix
        uniform mat4   u_view;          // View matrix
        uniform mat4   u_projection;    // Projection matrix
        attribute vec3 a_position;      // Vertex position
        void main()
        {
            // gl_Position = u_projection * u_view * u_model * vec4(a_position,1.0);
            gl_position = <transform>;
        }
        """

        fragment = """
        void main()
        {
            gl_FragColor = vec4(1.0, 0.0, 0.0, 1.0);
        }
        """

        #  vertex = """
        #  uniform vec4 u_color;
        #  attribute vec3 position;
        #  attribute vec4 color;
        #  varying vec4 v_color;
        #  void main()
        #  {
            #  v_color = u_color * color;
            #  gl_Position = <transform>;
        #  }
        #  """

        #  fragment = """
        #  varying vec4 v_color;
        #  void main()
        #  {
            #  gl_FragColor = v_color;
        #  }
        #  """

        window = app.Window(width=1024, height=1024,
                            color=(0.30, 0.30, 0.35, 1.00))


        verts = np.zeros(len(points), [("a_position", np.float32, 3)])
        verts["a_position"] = points


        verts = verts.view(gloo.VertexBuffer)
        indeces = np.array(triangle).astype(np.uint32)
        indeces = indeces.view(gloo.IndexBuffer)

        surface = gloo.Program(vertex, fragment)
        surface.bind(verts)
        surface['u_model'] = np.eye(4, dtype=np.float32)
        surface['u_view'] = glm.translation(0, 0, -5)
        surface['transform'] = Trackball(Position("a_position"))
        window.attach(surface['transform'])

        @window.event
        def on_draw(draw_triangles):
            """ draws and rotates the surface """
            global PHI, THETA
            window.clear()

            #  surface.draw(gl.GL_TRIANGLES, indeces)
            surface.draw(gl.GL_LINES, indeces)

            #  Make surface rotate
            THETA += 0.5 # degrees
            PHI += 0.5 # degrees
            model = np.eye(4, dtype=np.float32)
            glm.rotate(model, THETA, 0, 0, 1)
            glm.rotate(model, PHI, 0, 1, 0)
            surface['u_model'] = model


        @window.event
        def on_resize(width, height):
            """ deals with resizing the image """
            surface['u_projection'] = glm.perspective(45.0, width / float(height), 2.0, 100.0)

        @window.event
        def on_init():
            """ initializes opengl? """
            gl.glEnable(gl.GL_DEPTH_TEST)


        #  app.run(interactive=True)
        app.run()

def plot(data=None):
    """ simply calls the plot method """

    surface = GlumpyPlotter(data)
    surface.plot()
    
