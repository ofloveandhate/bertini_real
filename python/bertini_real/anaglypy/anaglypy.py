"""
Foong Min Wong
University of Wisconsin, Eau Claire
Summer 2019

.. module:: anaglypy
    :platform: Unix, Windows
    :synopsis: The anagylpy module contains python and shell scripts to generate 3d anaglyph video animation of algebraic surfaces.

"""
from os.path import dirname
import os
import bpy
import bertini_real
import json
import math
import numpy as np
import sys

# Delete any meshes before creating a new mesh
item = 'MESH'
bpy.ops.object.select_all(action='DESELECT')
bpy.ops.object.select_by_type(type=item)
bpy.ops.object.delete()
bpy.context.scene.render.use_overwrite = False

# Retrieve filename
fileName = os.getcwd().split(os.sep)[-1]

# Switch render mode to Cycle Render
bpy.context.scene.render.engine = 'CYCLES'

# Import movie data from JSON file
with open(dirname(os.getcwd())+'/' + sys.argv[-1]) as movie_data:
    prefs = json.load(movie_data)


class ReversableList(list):
    """ Create a ReversableList object for reversing order of data

        :param list: The list to be read.

    """

    def reverse(self):
        """ Reverse function for raw surface data

            :rtype: A reversed list

        """
        return list(reversed(self))


def extract_points(self):
    """ Helper method for plot_surface_samples()
        Extract points from vertices

        :param data: Surface decomposition data
        :rtype: List of tuples of length 3.
    """
    points = []

    for v in self.decomposition.vertices:
        q = [None] * 3

        for i in range(3):
            q[i] = v[i].real
        points.append(q)

    return points


# Define objects colors
r, g, b, alpha = prefs['color']['r'], prefs['color'][
    'g'], prefs['color']['b'], prefs['color']['alpha']


def diffuse():
    """ Set the general color of an object """

    # Set new material to variable
    mat = bpy.data.materials.new('MaterialName')

    bpy.data.materials['MaterialName'].use_nodes = True

    met = mat.node_tree.nodes["Principled BSDF"]

    nodes = mat.node_tree.nodes

    met.inputs[0].default_value = (r, g, b, alpha)

    met.inputs[4].default_value = 1

    node_principle = nodes.get("Principled BSDF")

    node_text = nodes.new(type='ShaderNodeTexCoord')

    mat.node_tree.links.new(
        node_text.outputs['Normal'], node_principle.inputs['Metallic'])

    # assign material to object
    bpy.context.object.data.materials.append(mat)


def diffuse2():

    # Set new material to variable
    mat = bpy.data.materials.new('MaterialName2')

    bpy.data.materials['MaterialName2'].use_nodes = True

    met = mat.node_tree.nodes["Principled BSDF"]

    nodes = mat.node_tree.nodes

    met.inputs[0].default_value = (r, g, b, alpha)

    met.inputs[4].default_value = 1

    node_principle = nodes.get("Principled BSDF")

    node_text = nodes.new(type='ShaderNodeTexCoord')

    mat.node_tree.links.new(
        node_text.outputs['Normal'], node_principle.inputs['Metallic'])

    # assign material to object
    bpy.context.object.data.materials.append(mat)


class Anaglypy():

    def __init__(self, data=None):
        """ Read data from disk

           :param data: Surface decomposition data. If data is None, then it reads the most recent BRData.pkl.
        """

        if data is None:
            self.decomposition = bertini_real.data.ReadMostRecent()
        else:
            self.decomposition = data

    def generate_fv(self):
        """ Generate faces and vertices of smooth surface """

        # Extract points
        points = extract_points(self)

        # Extract face
        face = self.decomposition.surface.surface_sampler_data

        # Define vertices, faces
        vertex = [p for p in points]
        faces = [y for x in face for y in x]

        return vertex, faces

    def generate_fv_raw(self):
        """ Generate faces and vertices of raw surface """

        points = extract_points(self)

        surf = self.decomposition.surface

        num_faces = surf.num_faces

        which_faces = list(range(num_faces))

        if not len(which_faces):
            which_faces = list(range(num_faces))

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

            T = []

            while 1:
                ## top edge ##
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

                    if(curr_edge == -10):
                        continue

                    if (curr_edge[0] < 0 and curr_edge[1] < 0 and curr_edge[2] < 0):
                        continue

                    curr_edge = ReversableList(curr_edge)
                    curr_edge = curr_edge.reverse()

                ## bottom edge ##
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

                    if(curr_edge == -10):
                        continue

                    if (curr_edge[0] < 0 and curr_edge[1] < 0 and curr_edge[2] < 0):
                        continue

                ## left edge ##
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

                ## right edge ##
                elif case == 4:

                    if right_edge_counter < face['num right']:

                        if face['right'][right_edge_counter] < 0:
                            continue

                        slice_ind = face['middle slice index'] + 1
                        edge_ind = face['right'][right_edge_counter]
                        curr_edge = surf.critical_point_slices[
                            slice_ind].edges[edge_ind]
                        right_edge_counter = right_edge_counter + 1

                        curr_edge = ReversableList(curr_edge)
                        curr_edge = curr_edge.reverse()

                    else:
                        case += 1
                        continue

                ## last case ##
                elif case == 5:
                    break

                t1 = [points[curr_edge[0]], points[curr_edge[1]],
                      points[face['midpoint']]]
                t2 = [points[curr_edge[1]], points[curr_edge[2]],
                      points[face['midpoint']]]

                t3 = (curr_edge[0], curr_edge[1], face['midpoint'])
                t4 = (curr_edge[1], curr_edge[2], face['midpoint'])

                T.append(t1)
                T.append(t2)

                TT.append(t3)
                TT.append(t4)

        # Define vertices, faces
        vertex = [p for p in points]
        faces = TT

        return vertex, faces

    def generate_obj_scene(self, vertex, faces):
        """ Generate scene for one mesh
            :param vertex: vertices of a mesh
            :param faces: faces of a mesh
        """

        # Define mesh and object's name
        mesh = bpy.data.meshes.new(fileName)
        object = bpy.data.objects.new(fileName, mesh)

        # Set location and scene of object
        object.location = bpy.context.scene.cursor.location
        bpy.context.scene.collection.objects.link(object)

        # Create mesh (error here, it works)
        mesh.from_pydata(vertex, [], faces)
        mesh.update()

        # Make object active
        bpy.context.view_layer.objects.active = object

        diffuse()

        # Retrieve object dimensions
        object_dimensions = object.dimensions

        # Resize/ Scale object
        bpy.context.object.dimensions = object.dimensions[0], prefs['object']['dimension'], object.dimensions[
            2]  # resize z to 3

        # Grab the current object scale
        object_scale = object.scale

        # Scale them by z scale
        object.scale = (object_scale[1], object_scale[1], object_scale[1])

        # Rescale them (should try ratio method?)
        object.scale = (object_scale[1] + prefs['object']['inflation'],
                        object_scale[1] + prefs['object']['inflation'], object_scale[1] + prefs['object']['inflation'])

        # go edit mode
        bpy.ops.object.mode_set(mode='EDIT')

        # select all faces
        bpy.ops.mesh.select_all(action='SELECT')

        # recalculate outside normals
        bpy.ops.mesh.normals_make_consistent(inside=False)

        # go object mode again
        bpy.ops.object.editmode_toggle()

        context = bpy.context
        scene = context.scene

        return object, scene

    def generate_both_scene(self, vertex, faces, vertex_raw, faces_raw):
        """ Generate scene for two objects

            :param vertex: vertices of a smooth mesh
            :param faces: faces of a smooth mesh
            :param vertex_raw: vertices of a raw mesh
            :param faces_raw: faces of a raw mesh
        """

        # Define mesh and object's name
        mesh = bpy.data.meshes.new(fileName)
        object = bpy.data.objects.new(fileName, mesh)

        # Set location and scene of object
        object.location = bpy.context.scene.cursor.location
        object.location = (prefs['object_both']['location'], prefs[
                           'object_both']['location'], 0)
        bpy.context.scene.collection.objects.link(object)

        # Create mesh
        mesh.from_pydata(vertex, [], faces)
        mesh.update()

        # Make object active
        bpy.context.view_layer.objects.active = object

        diffuse()

        # Retrieve object dimensions
        object_dimensions = object.dimensions

        # Resize/ Scale object
        bpy.context.object.dimensions = object.dimensions[0], prefs['object_both']['dimension'], object.dimensions[
            2]

        # Grab the current object scale
        object_scale = object.scale

        # Scale them by z scale
        object.scale = (object_scale[1], object_scale[1], object_scale[1])

        # Rescale them (should try ratio method?)
        object.scale = (object_scale[1] + prefs['object_both']['inflation'],
                        object_scale[1] + prefs['object_both']['inflation'], object_scale[1] + prefs['object_both']['inflation'])

        # go edit mode
        bpy.ops.object.mode_set(mode='EDIT')

        # select all faces
        bpy.ops.mesh.select_all(action='SELECT')

        # recalculate outside normals
        bpy.ops.mesh.normals_make_consistent(inside=False)

        # Define 2nd mesh and object's name
        mesh1 = bpy.data.meshes.new(fileName)
        object1 = bpy.data.objects.new(fileName, mesh1)

        # Set location and scene of object
        object1.location = (-prefs['object_both']
                            ['location'], -prefs['object_both']['location'], 0)
        bpy.context.scene.collection.objects.link(object1)

        # Create mesh
        mesh1.from_pydata(vertex_raw, [], faces_raw)
        mesh1.update()

        # Make object active
        bpy.context.view_layer.objects.active = object1

        diffuse2()

        # Retrieve object dimensions
        object1_dimensions = object1.dimensions

        # Resize/ Scale object
        bpy.context.object.dimensions = object1.dimensions[0], prefs['object_both']['dimension'], object1.dimensions[
            2]

        # Grab the current object scale
        object1_scale = object1.scale

        # Scale them by z scale
        object1.scale = (object1_scale[1], object1_scale[1], object1_scale[1])

        # Rescale them (should try ratio method?)
        object1.scale = (object1_scale[1] + prefs['object_both']['inflation'],
                         object1_scale[1] + prefs['object_both']['inflation'], object1_scale[1] + prefs['object_both']['inflation'])

        # go edit mode
        bpy.ops.object.mode_set(mode='EDIT')

        # select all faces
        bpy.ops.mesh.select_all(action='SELECT')

        # recalculate outside normals
        bpy.ops.mesh.normals_make_consistent(inside=False)

        context = bpy.context
        scene = context.scene

        return object, object1, scene

    def generate_multi_obj_scene(self, vertex, faces):
        """ Generate scene for three objects 

            :param vertex: vertices of a smooth mesh
            :param faces: faces of a smooth mesh

        """

        # Define mesh and object's name
        mesh = bpy.data.meshes.new(fileName)
        object = bpy.data.objects.new(fileName, mesh)

        # Set location and scene of object
        object.location = bpy.context.scene.cursor.location
        bpy.context.scene.collection.objects.link(object)

        # Create mesh (error here, it works)
        mesh.from_pydata(vertex, [], faces)
        mesh.update()

        # Make object active
        bpy.context.view_layer.objects.active = object

        diffuse()

        # Retrieve object dimensions
        object_dimensions = object.dimensions

        # Resize/ Scale object
        bpy.context.object.dimensions = object.dimensions[
            0], prefs['object_multi']['dimension'], object.dimensions[2]

        # Grab the current object scale
        object_scale = object.scale

        # Scale them by z scale
        object.scale = (object_scale[1], object_scale[1], object_scale[1])

        # Rescale them (should try ratio method?)
        object.scale = (object_scale[1] + prefs['object_multi']['inflation'],
                        object_scale[1] + prefs['object_multi']['inflation'], object_scale[1] + prefs['object_multi']['inflation'])

        # go edit mode
        bpy.ops.object.mode_set(mode='EDIT')

        # select all faces
        bpy.ops.mesh.select_all(action='SELECT')

        # recalculate outside normals
        bpy.ops.mesh.normals_make_consistent(inside=False)

        object1 = bpy.data.objects.new(object.name, object.data.copy())

        # adds the object to the active scene
        bpy.context.collection.objects.link(object1)

        # Make object active
        bpy.context.view_layer.objects.active = object1

        # Retrieve object dimensions
        object1_dimensions = object1.dimensions

        # Resize/ Scale object
        bpy.context.object.dimensions = object.dimensions[0], prefs['object_multi']['dimension'], object.dimensions[
            2]

        # Grab the current object scale
        object1_scale = object1.scale

        # Scale them by z scale
        object1.scale = (object1_scale[1], object1_scale[1], object1_scale[1])

        # Rescale them (should try ratio method?)
        object1.scale = (object1_scale[1] + prefs['object_multi']['inflation'],
                         object1_scale[1] + prefs['object_multi']['inflation'], object1_scale[1] + prefs['object_multi']['inflation'])

        object1.location = (prefs['object_multi']['location'], prefs[
                            'object_multi']['location'], 0)

        object2 = bpy.data.objects.new(object.name, object.data.copy())

        # adds the object to the active scene
        bpy.context.collection.objects.link(object2)

        # Make object active
        bpy.context.view_layer.objects.active = object2

        # Retrieve object dimensions
        object2_dimensions = object2.dimensions

        # Resize/ Scale object
        bpy.context.object.dimensions = object.dimensions[0], prefs['object_multi']['dimension'], object.dimensions[
            2]  # resize z to 1.0

        # Grab the current object scale
        object2_scale = object2.scale

        # Scale them by z scale
        object2.scale = (object2_scale[1], object2_scale[1], object2_scale[1])

        # Rescale them (should try ratio method?)
        object2.scale = (object2_scale[1] + prefs['object_multi']['inflation'],
                         object2_scale[1] + prefs['object_multi']['inflation'], object2_scale[1] + prefs['object_multi']['inflation'])

        object2.location = (-prefs['object_multi']
                            ['location'], -prefs['object_multi']['location'], 0)

        # go object mode again
        bpy.ops.object.editmode_toggle()

        context = bpy.context
        scene = context.scene

        # bpy.context.scene.render.image_settings.views_format = 'STEREO_3D'

        # bpy.data.cameras['Camera'].stereo.convergence_distance = 11
        # bpy.data.cameras['Camera'].stereo.interocular_distance = 1.5

        return object, object1, object2, scene

    def rotate_z(self, object, scene):
        """ Create a rotation around the object z-axis

            :param object: meshes of faces, edges and/or vertices
            :param scene: place to store objects and materials
        """

        scene.frame_start = 0
        scene.frame_end = prefs['video_settings']['num_frames']

        # rotate nothing
        object.rotation_euler = (0.0, 0.0, 0.0)
        object.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the z-axis
        object.rotation_euler = (0, 0, math.pi)
        object.keyframe_insert(data_path='rotation_euler',
                               frame=scene.frame_end // 2)

        # rotate at the z-axis
        object.rotation_euler = (0, 0, math.pi * 2)
        object.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end)

    def rotate_xyz(self, object, scene):
        """ Create a rotation around the object x, y and, z-axis

            :param object: meshes of faces, edges and/or vertices
            :param scene: place to store objects and materials
        """

        scene.frame_start = 0
        scene.frame_end = prefs['video_settings']['num_frames']

        # rotate nothing
        object.rotation_euler = (0.0, 0.0, 0.0)
        object.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the z-axis
        object.rotation_euler = (0, 0, math.pi * 2)
        object.keyframe_insert(data_path='rotation_euler',
                               frame=scene.frame_end // 3)

        # rotate at the y-axis
        object.rotation_euler = (0, math.pi * 2, math.pi * 2)
        object.keyframe_insert(data_path='rotation_euler',
                               frame=scene.frame_end // 1.5)

        # rotate at the x-axis
        object.rotation_euler = (math.pi * 2, math.pi * 2, math.pi * 2)
        object.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end)

    def spin(self, object, scene):
        """ Create a spinning animation of object 

            :param object: meshes of faces, edges and/or vertices
            :param scene: place to store objects and materials

        """

        scene.frame_start = 0
        scene.frame_end = prefs['video_settings']['num_frames']

        # rotate nothing
        object.rotation_euler = (0.0, 0.0, 0.0)
        object.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the z-axis
        object.rotation_euler = (0, 0, math.pi * 2)
        object.keyframe_insert(data_path='rotation_euler',
                               frame=scene.frame_end // 2)

        # rotate at the z-axis
        object.rotation_euler = (0, 0, 0)
        object.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end)

    def rotate_z_both(self, object, object1, scene):
        """ Create a rotation around two objects z-axis 

            :param object: first meshes of faces, edges and/or vertices
            :param object1: second meshes of faces, edges and/or vertices
            :param scene: place to store objects and materials
        """

        scene.frame_start = 0
        scene.frame_end = prefs['video_settings']['num_frames']

        # rotate nothing
        object.rotation_euler = (0.0, 0.0, 0.0)
        object.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the z-axis
        object.rotation_euler = (0, 0, math.pi)
        object.keyframe_insert(data_path='rotation_euler',
                               frame=scene.frame_end // 2)

        # rotate at the z-axis
        object.rotation_euler = (0, 0, math.pi * 2)
        object.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end)

        # rotate nothing
        object1.rotation_euler = (0.0, 0.0, 0.0)
        object1.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the z-axis
        object1.rotation_euler = (0, 0, math.pi)
        object1.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end // 2)

        # rotate at the z-axis
        object1.rotation_euler = (0, 0, math.pi * 2)
        object1.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end)

    def rotate_xyz_both(self, object, object1, scene):
        """ Create a rotation around two objects x, y and, z-axis 

            :param object: first meshes of faces, edges and/or vertices
            :param object1: second meshes of faces, edges and/or vertices
            :param scene: place to store objects and materials

        """

        scene.frame_start = 0
        scene.frame_end = prefs['video_settings']['num_frames']

        # rotate nothing
        object.rotation_euler = (0.0, 0.0, 0.0)
        object.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the z-axis
        object.rotation_euler = (0, 0, math.pi * 2)
        object.keyframe_insert(data_path='rotation_euler',
                               frame=scene.frame_end // 3)

        # rotate at the y-axis
        object.rotation_euler = (0, math.pi * 2, math.pi * 2)
        object.keyframe_insert(data_path='rotation_euler',
                               frame=scene.frame_end // 1.5)

        # rotate at the x-axis
        object.rotation_euler = (math.pi * 2, math.pi * 2, math.pi * 2)
        object.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end)

        object1.rotation_euler = (0.0, 0.0, 0.0)
        object1.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the z-axis
        object1.rotation_euler = (0, 0, math.pi * 2)
        object1.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end // 3)

        # rotate at the y-axis
        object1.rotation_euler = (0, math.pi * 2, math.pi * 2)
        object1.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end // 1.5)

        # rotate at the x-axis
        object1.rotation_euler = (math.pi * 2, math.pi * 2, math.pi * 2)
        object1.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end)

    def spin_both(self, object, object1, scene):
        """ Create a spinning animation of two objects

            :param object: first meshes of faces, edges and/or vertices
            :param object1: second meshes of faces, edges and/or vertices
            :param scene: place to store objects and materials

        """

        scene.frame_start = 0
        scene.frame_end = prefs['video_settings']['num_frames']

        # rotate nothing
        object.rotation_euler = (0.0, 0.0, 0.0)
        object.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the z-axis
        object.rotation_euler = (0, 0, math.pi * 2)
        object.keyframe_insert(data_path='rotation_euler',
                               frame=scene.frame_end // 2)

        # rotate at the z-axis
        object.rotation_euler = (0, 0, 0)
        object.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end)

        # rotate nothing
        object1.rotation_euler = (0.0, 0.0, 0.0)
        object1.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the z-axis
        object1.rotation_euler = (0, 0, math.pi * 2)
        object1.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end // 2)

        # rotate at the z-axis
        object1.rotation_euler = (0, 0, 0)
        object1.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end)

    def multi_rotate(self, object, object1, object2, scene):
        """ Create XYZ rotation animation for three objects 

            :param object: first meshes of faces, edges and/or vertices
            :param object1: second meshes of faces, edges and/or vertices
            :param object2: third meshes of faces, edges and/or vertices
            :param scene: place to store objects and materials

        """

        scene.frame_start = 0
        scene.frame_end = prefs['video_settings']['num_frames']

        # rotate nothing
        object.rotation_euler = (0.0, 0.0, 0.0)
        object.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the y-axis
        object.rotation_euler = (0, math.pi, 0)
        object.keyframe_insert(data_path='rotation_euler',
                               frame=scene.frame_end // 2)

        # rotate at the y-axis
        object.rotation_euler = (0, math.pi * 2, 0)
        object.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end)

        # rotate nothing
        object1.rotation_euler = (0.0, 0.0, 0.0)
        object1.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the z-axis
        object1.rotation_euler = (0, 0, math.pi)
        object1.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end // 2)

        # rotate at the z-axis
        object1.rotation_euler = (0, 0, math.pi * 2)
        object1.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end)

        # rotate nothing
        object2.rotation_euler = (0.0, 0.0, 0.0)
        object2.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the x-axis
        object2.rotation_euler = (math.pi, 0, 0)
        object2.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end // 2)

        # rotate at the x-axis
        object2.rotation_euler = (math.pi * 2, 0, 0)
        object2.keyframe_insert(
            data_path='rotation_euler', frame=scene.frame_end)


def stereo_flag(scene, directory):
    """ Create Anaglyph 3D or Non-Anaglyph videos

        :param scene: place to store objects and materials
        :param directory: directory of surface
    """

    if prefs['anaglyph'] == 1:
        scene.render.use_multiview = True
        bpy.context.scene.render.image_settings.views_format = prefs[
            'video_settings']['views_format']

        bpy.data.cameras['Camera'].stereo.convergence_distance = prefs[
            'convergence_distance']
        bpy.data.cameras['Camera'].stereo.interocular_distance = prefs[
            'interocular_distance']

        directory = directory + "_anaglyph"
    else:
        scene.render.use_multiview = False
        directory = directory + "_non_anaglyph"

    return directory


def render(scene, directory):
    """ Create XYZ rotation animation for three objects 

        :param scene: place to store objects and materials
        :param directory: directory of surface

    """

    scene.render.resolution_percentage = prefs['resolution']['percentage']
    scene.render.resolution_x = prefs['resolution']['x']
    scene.render.resolution_y = prefs['resolution']['y']

    directory = stereo_flag(scene, directory)

    scene.render.filepath = "render/" + fileName + directory
    scene.render.image_settings.file_format = prefs[
        'video_settings']['file_format']

    scene.render.image_settings.color_mode = prefs[
        'video_settings']['color_mode']

    bpy.context.scene.cycles.film_exposure = prefs[
        'video_settings']['film_exposure']

    bpy.ops.render.render(animation=True)

    print("Export " + '\x1b[0;33;40m' + "Anaglyph 3D " + '\x1b[0m' +
          '\x1b[0;35;40m' + fileName + directory + ".avi" + '\x1b[0m' + " successfully")

    bpy.ops.wm.save_as_mainfile(
        filepath=os.getcwd() + "/render/" + fileName + directory + ".blend")


def raw(data=None):
    """ Generate animation of raw surface """
    surface = Anaglypy(data)
    vertex, faces = surface.generate_fv_raw()

    if prefs['animation'] == "Rotate_Z":
        object, scene = surface.generate_obj_scene(vertex, faces)
        directory = "_rotate_z_raw"
        surface.rotate_z(object, scene)
        render(scene, directory)

    elif prefs['animation'] == "Rotate_XYZ":
        object, scene = surface.generate_obj_scene(vertex, faces)
        directory = "_rotate_xyz_raw"
        surface.rotate_xyz(object, scene)
        render(scene, directory)

    elif prefs['animation'] == "Spin":
        object, scene = surface.generate_obj_scene(vertex, faces)
        directory = "_spin_raw"
        surface.spin(object, scene)
        render(scene, directory)

    elif prefs['animation'] == "Multi_Rotate":
        directory = "_multi_rotate_raw"
        object, object1, object2, scene = surface.generate_multi_obj_scene(
            vertex, faces)
        surface.multi_rotate(object, object1, object2, scene)
        render(scene, directory)


def smooth(data=None):
    """ Generate animation of smooth surface """

    surface = Anaglypy(data)
    vertex, faces = surface.generate_fv()

    if prefs['animation'] == "Rotate_Z":
        object, scene = surface.generate_obj_scene(vertex, faces)
        directory = "_rotate_z_smooth"
        surface.rotate_z(object, scene)
        render(scene, directory)

    elif prefs['animation'] == "Rotate_XYZ":
        object, scene = surface.generate_obj_scene(vertex, faces)
        directory = "_rotate_xyz_smooth"
        surface.rotate_xyz(object, scene)
        render(scene, directory)

    elif prefs['animation'] == "Spin":
        object, scene = surface.generate_obj_scene(vertex, faces)
        directory = "_spin_smooth"
        surface.spin(object, scene)
        render(scene, directory)

    elif prefs['animation'] == "Multi_Rotate":
        object, object1, object2, scene = surface.generate_multi_obj_scene(
            vertex, faces)
        directory = "_multi_rotate_smooth"
        surface.multi_rotate(object, object1, object2, scene)
        render(scene, directory)


def both(data=None):
    """ Generate animation of raw and smooth surfaces """

    surface = Anaglypy(data)
    vertex, faces = surface.generate_fv()
    vertex_raw, faces_raw = surface.generate_fv_raw()

    if prefs['animation'] == "Rotate_Z":
        object, object1, scene = surface.generate_both_scene(
            vertex, faces, vertex_raw, faces_raw)
        directory = "_rotate_z_both"
        surface.rotate_z_both(object, object1, scene)
        render(scene, directory)

    elif prefs['animation'] == "Rotate_XYZ":
        object, object1, scene = surface.generate_both_scene(
            vertex, faces, vertex_raw, faces_raw)
        directory = "_rotate_xyz_both"
        surface.rotate_xyz_both(object, object1, scene)
        render(scene, directory)

    elif prefs['animation'] == "Spin":
        object, object1, scene = surface.generate_both_scene(
            vertex, faces, vertex_raw, faces_raw)
        directory = "_spin_both"
        surface.spin_both(object, object1, scene)
        render(scene, directory)


if __name__ == "__main__":
    if prefs['raw_smooth'] == "Raw":
        raw()
    elif prefs['raw_smooth'] == "Smooth":
        smooth()
    elif prefs['raw_smooth'] == "Both":
        both()

bpy.ops.wm.quit_blender()
bpy.ops.wm.window_close()
