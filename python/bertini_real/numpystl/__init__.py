"""
Foong Min Wong
University of Wisconsin, Eau Claire
Fall 2018 - Spring 2019
Implementing stereolithography (STL) export feature for Bertini_real

Current Version:


TODO:

"""
import os
import numpy as np
from bertini_real.data import BRData
from bertini_real.surface import Surface, Curve
import bertini_real
import bertini_real.util
import dill
import numpy as np
import matplotlib
from stl import mesh
import trimesh


class NumpySTL():
    """ creates the glumpyplotter object """

    def __init__(self, data=None):
        if data is None:
            self.decomposition = bertini_real.data.ReadMostRecent()
        else:
            self.decomposition = data

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

    def stl_raw(self):
    	print("stl_raw in progress")
    	points = self.extract_points()
    	surf = self.decomposition.surface



    def stl_sampler(self):
        points = self.extract_points()
        faces = self.decomposition.surface.surface_sampler_data
        # add vertex and surface to mesh
        vertex = []

        for p in points:
            vertex.append(p)

        vertex_np_array = np.array(vertex)
        face = []

        for f in faces:  # for each face
            for tri in f:  # for triangle in face
                face.append([tri[0], tri[1], tri[2]])

        face_np_array = np.array(face)

        obj = mesh.Mesh(np.zeros(face_np_array.shape[
                        0], dtype=mesh.Mesh.dtype))

        for i, f in enumerate(face_np_array):
            for j in range(3):
                obj.vectors[i][j] = vertex_np_array[f[j], :]

        # get object filename
        fileName = os.getcwd().split(os.sep)[-1]

        obj.save('a' + fileName + '.stl')

        normmesh = trimesh.load_mesh('a' + fileName + '.stl')

        normmesh.fix_normals()

        for facet in normmesh.facets:
            normmesh.visual.face_colors[facet] = trimesh.visual.random_color()

        # normmesh.show()

        

        normmesh.export(file_obj='anorm' + fileName + '.stl', file_type='stl')

        print("Export " + '\x1b[2;30;45m' + "anorm" + fileName + ".stl" + '\x1b[0m' + " successfully")


# ------------------------------------------------------------------------------------- #


# ------------------------------------------------------------------------------------- #

def stl_sampler(data=None):
    surface = NumpySTL(data)
    surface.stl_sampler()

def stl_raw(data=None):
    surface = NumpySTL(data)
    surface.stl_raw()
