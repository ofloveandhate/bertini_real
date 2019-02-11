"""
Dan Hessler
University of Wisconsin, Eau Claire
Fall 2018 - Spring 2019
Porting to Glumpy for faster surface rendering
using OpenGL

Foong Min Wong
University of Wisconsin, Eau Claire
Fall 2018 - Spring 2019
Implementing stereolithography export feature

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

import os
import numpy as np
from glumpy import app, gl, glm, gloo
from glumpy.transforms import Trackball, Position
import bertini_real
from stl import mesh
import trimesh

class GlumpyPlotter():
    """ creates the glumpyplotter object """
    def __init__(self, data=None):
        if data is None:
            self.decomposition = bertini_real.data.ReadMostRecent()
        else:
            self.decomposition = data

    def fvtostl(self):

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

        data = self.decomposition
        points = extract_points(data)
        faces = self.decomposition.surface.surface_sampler_data
       
        # add vertex and surface to mesh
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

        # get object filename
        fileName = os.getcwd().split(os.sep)[-1]

        obj.save('a' + fileName + '.stl')

        normmesh = trimesh.load_mesh('a' + fileName + '.stl')

        normmesh.fix_normals()

        for facet in normmesh.facets:
            normmesh.visual.face_colors[facet] = trimesh.visual.random_color()

        # normmesh.show()
        
        normmesh.export(file_obj='anorm' + fileName + '.stl', file_type='stl')
        print('Export STL successfully!')

    def plot(self, color_function):
        """ method used to plot a surface """
        print("Plotting object of dimension: {}".format(self.decomposition.dimension))

        data = self.decomposition
        faces = data.surface.surface_sampler_data

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

        def make_colors(points):
            """
            computes colors according to a function
            """

            colors = []

            for i in range(len(points)):
                x = points[i][0]
                y = points[i][1]
                z = points[i][2]

                f1, f2, f3 = color_function

                r = f1(x,y,z)
                g = f2(x,y,z)
                b = f3(x,y,z)
                colors.append([r, g, b, 1])
            return colors

        points = extract_points(data)

        triangle = []

        for i in range(len(faces)):
            for tri in faces[i]:
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

def plot(color_function=None, data=None):
    """
    simply calls the plot method
    color_function contains 3 functions to compute the colors
    of the surface
    """
    # default colors if none are provided
    # function = x^2 + y^2 + z^2
    # sort of
    if color_function == None:
        def f1(x,y,z):
            return x**2
        def f2(x,y,z):
            return y**2
        def f3(x,y,z):
            return z**2

        color_function = f1, f2, f3

    surface = GlumpyPlotter(data)
    surface.plot(color_function)

def fvtostl(data=None):
    surface = GlumpyPlotter(data)
    surface.fvtostl()
