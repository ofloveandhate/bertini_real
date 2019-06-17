import os
import bpy
import bertini_real
import math
import numpy as np

# Delete any meshes before creating a new mesh
item = 'MESH'
bpy.ops.object.select_all(action='DESELECT')
bpy.ops.object.select_by_type(type=item)
bpy.ops.object.delete()
bpy.context.scene.render.use_overwrite = False

# Retrieve filename
fileName = os.getcwd().split(os.sep)[-1]

bpy.context.scene.render.engine = 'CYCLES'


class ReversableList(list):
    """ Create a ReversableList object for reversing order of data

        :param list: The list to be read.

    """

    def reverse(self):
        """ Reverse function for raw surface data

            :rtype: A reversed list

        """
        return list(reversed(self))


def user_pick(options):
    """ Create a list of user choice options

        :param options: The list of options

    """

    print("Please choose:")

    for idx, element in enumerate(options):
        print("{}) {}".format(idx + 1, element))
    i = input("Enter number: ")
    print("\n")
    try:
        if 0 < int(i) <= len(options):
            return int(i)
    except:
        pass
    return None


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
            q[i] = v['point'][i].real
        points.append(q)

    return points


class Anaglypy():

    def __init__(self, data=None):

        if data is None:
            self.decomposition = bertini_real.data.ReadMostRecent()
        else:
            self.decomposition = data

    def generate_fv(self):

        # Extract points
        points = extract_points(self)

        # Extract face
        face = self.decomposition.surface.surface_sampler_data

        # Define vertices, faces
        vertex = [p for p in points]
        faces = [y for x in face for y in x]

        return vertex, faces

    def generate_fv_raw(self):
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

        # Retrieve object dimensions
        object_dimensions = object.dimensions

        # Resize/ Scale object
        bpy.context.object.dimensions = object.dimensions[
            0], object.dimensions[1], 1.5  # resize z to 1.5

        # Grab the current object scale
        object_scale = object.scale

        # Scale them by z scale
        object.scale = (object_scale[2], object_scale[2], object_scale[2])

        # Rescale them (should try ratio method?)
        object.scale = (object_scale[2] + 0.2,
                        object_scale[2] + 0.2, object_scale[2] + 0.2)

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

        scene.render.use_multiview = True

        bpy.context.scene.render.image_settings.views_format = 'STEREO_3D'

        bpy.context.scene.cycles.film_exposure = 4.5

        bpy.data.cameras['Camera'].stereo.convergence_distance = 15
        bpy.data.cameras['Camera'].stereo.interocular_distance = 0.15

        return object, scene

    def generate_multi_obj_scene(self, vertex, faces):

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

        # Retrieve object dimensions
        object_dimensions = object.dimensions

        # Resize/ Scale object
        bpy.context.object.dimensions = object.dimensions[
            0], object.dimensions[1], 1.0  # resize z to 1.0

        # Grab the current object scale
        object_scale = object.scale

        # Scale them by z scale
        object.scale = (object_scale[2], object_scale[2], object_scale[2])

        # Rescale them (should try ratio method?)
        object.scale = (object_scale[2] + 0.15,
                        object_scale[2] + 0.15, object_scale[2] + 0.15)

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
        bpy.context.object.dimensions = object.dimensions[
            0], object.dimensions[1], 1.0  # resize z to 1.0

        # Grab the current object scale
        object1_scale = object1.scale

        # Scale them by z scale
        object1.scale = (object1_scale[2], object1_scale[2], object1_scale[2])

        # Rescale them (should try ratio method?)
        object1.scale = (object1_scale[2] + 0.15,
                         object1_scale[2] + 0.15, object1_scale[2] + 0.15)

        object1.location = (1.5, 1.5, 0)

        object2 = bpy.data.objects.new(object.name, object.data.copy())

        # adds the object to the active scene
        bpy.context.collection.objects.link(object2)

        # Make object active
        bpy.context.view_layer.objects.active = object2

        # Retrieve object dimensions
        object2_dimensions = object2.dimensions

        # Resize/ Scale object
        bpy.context.object.dimensions = object.dimensions[
            0], object.dimensions[1], 1.0  # resize z to 1.0

        # Grab the current object scale
        object2_scale = object2.scale

        # Scale them by z scale
        object2.scale = (object2_scale[2], object2_scale[2], object2_scale[2])

        # Rescale them (should try ratio method?)
        object2.scale = (object2_scale[2] + 0.15,
                         object2_scale[2] + 0.15, object2_scale[2] + 0.15)

        object2.location = (-1.5, -1.5, 0)

        # go object mode again
        bpy.ops.object.editmode_toggle()

        context = bpy.context
        scene = context.scene

        scene.render.use_multiview = True

        bpy.context.scene.render.image_settings.views_format = 'STEREO_3D'

        bpy.context.scene.cycles.film_exposure = 4.5

        bpy.data.cameras['Camera'].stereo.convergence_distance = 15
        bpy.data.cameras['Camera'].stereo.interocular_distance = 0.15

        return object, object1, object2, scene

    def rotate_z(self, object, scene):
        # object.rotation_mode = 'XYZ'

        scene.frame_start = 0
        scene.frame_end = 200

        # rotate nothing
        object.rotation_euler = (0.0, 0.0, 0.0)
        object.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the z-axis
        object.rotation_euler = (0, 0, math.pi)
        object.keyframe_insert(data_path='rotation_euler', frame=100)

        # rotate at the z-axis
        object.rotation_euler = (0, 0, math.pi * 2)
        object.keyframe_insert(data_path='rotation_euler', frame=200)

    def rotate_xyz(self, object, scene):

        scene.frame_start = 0
        scene.frame_end = 300

        # rotate nothing
        object.rotation_euler = (0.0, 0.0, 0.0)
        object.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the z-axis
        object.rotation_euler = (0, 0, math.pi * 2)
        object.keyframe_insert(data_path='rotation_euler', frame=100)

        # rotate at the y-axis
        object.rotation_euler = (0, math.pi * 2, math.pi * 2)
        object.keyframe_insert(data_path='rotation_euler', frame=200)

        # rotate at the x-axis
        object.rotation_euler = (math.pi * 2, math.pi * 2, math.pi * 2)
        object.keyframe_insert(data_path='rotation_euler', frame=300)

    def spin_bf(self, object, scene):

        scene.frame_start = 0
        scene.frame_end = 200

        # rotate nothing
        object.rotation_euler = (0.0, 0.0, 0.0)
        object.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the z-axis
        object.rotation_euler = (0, 0, math.pi * 2)
        object.keyframe_insert(data_path='rotation_euler', frame=100)

        # rotate at the z-axis
        object.rotation_euler = (0, 0, 0)
        object.keyframe_insert(data_path='rotation_euler', frame=200)

    def multi_rotate(self, object, object1, object2, scene):

        scene.frame_start = 0
        scene.frame_end = 150

        # rotate nothing
        object.rotation_euler = (0.0, 0.0, 0.0)
        object.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the y-axis
        object.rotation_euler = (0, math.pi, 0)
        object.keyframe_insert(data_path='rotation_euler', frame=50)

        # rotate at the y-axis
        object.rotation_euler = (0, math.pi * 2, 0)
        object.keyframe_insert(data_path='rotation_euler', frame=150)

        # rotate nothing
        object1.rotation_euler = (0.0, 0.0, 0.0)
        object1.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the z-axis
        object1.rotation_euler = (0, 0, math.pi)
        object1.keyframe_insert(data_path='rotation_euler', frame=50)

        # rotate at the z-axis
        object1.rotation_euler = (0, 0, math.pi * 2)
        object1.keyframe_insert(data_path='rotation_euler', frame=150)

        # rotate nothing
        object2.rotation_euler = (0.0, 0.0, 0.0)
        object2.keyframe_insert(data_path='rotation_euler', frame=0)

        # rotate at the x-axis
        object2.rotation_euler = (math.pi, 0, 0)
        object2.keyframe_insert(data_path='rotation_euler', frame=50)

        # rotate at the x-axis
        object2.rotation_euler = (math.pi * 2, 0, 0)
        object2.keyframe_insert(data_path='rotation_euler', frame=150)

    def translate(self, object, scene):
        frame_num = 0

        positions = (0, 3, 2), (4, 1, 5), (3, -3, 1), (3, 3, 1), (1, 4, 1)
        start_pos = (0, 0, 0)

        for pos in positions:
            bpy.context.scene.frame_set(frame_num)
            object.location = pos
            object.keyframe_insert(data_path="location", index=-1)
            frame_num += 20


def render(scene, directory):
    scene.render.resolution_percentage = 100
    scene.render.resolution_x = 800
    scene.render.resolution_y = 600

    scene.render.filepath = "render/" + fileName + directory
    scene.render.image_settings.file_format = "AVI_JPEG"
    scene.render.image_settings.color_mode = "RGB"

    bpy.ops.render.render(animation=True)

    print("Export " + '\x1b[0;33;40m' + "Anaglyph 3D " + '\x1b[0m' +
          '\x1b[0;35;40m' + fileName + directory + ".avi" + '\x1b[0m' + " successfully")

    bpy.ops.wm.save_as_mainfile(
        filepath=os.getcwd() + "/render/" + fileName + directory + ".blend")


options = ["Rotate Z", "Rotate XYZ",
           "Spin Back & Forth", "Multi Rotate", "Translate"]


def smooth(data=None):
    surface = Anaglypy(data)
    vertex, faces = surface.generate_fv()

    option = user_pick(options)

    # wanna do switch/case or if/else statement?

    if option == 1:
        object, scene = surface.generate_obj_scene(vertex, faces)
        directory = "_rotate_z_smooth"
        surface.rotate_z(object, scene)
        render(scene, directory)

    elif option == 2:
        object, scene = surface.generate_obj_scene(vertex, faces)
        directory = "_rotate_xyz_smooth"
        surface.rotate_xyz(object, scene)
        render(scene, directory)

    elif option == 3:
        object, scene = surface.generate_obj_scene(vertex, faces)
        directory = "_spin_smooth"
        surface.spin_bf(object, scene)
        render(scene, directory)

    elif option == 4:
        object, object1, object2, scene = surface.generate_multi_obj_scene(
            vertex, faces)
        directory = "_multi_rotate_smooth"
        surface.multi_rotate(object, object1, object2, scene)
        render(scene, directory)

    elif option == 5:
        object, scene = surface.generate_obj_scene(vertex, faces)
        directory = "_rwalk_smooth"
        surface.translate(object, scene)
        render(scene, directory)


def raw(data=None):
    surface = Anaglypy(data)
    vertex, faces = surface.generate_fv_raw()

    option = user_pick(options)

    # wanna do switch/case or if/else statement?

    if option == 1:
        object, scene = surface.generate_obj_scene(vertex, faces)
        directory = "_rotate_z_raw"
        surface.rotate_z(object, scene)
        render(scene, directory)

    elif option == 2:
        object, scene = surface.generate_obj_scene(vertex, faces)
        directory = "_rotate_xyz_raw"
        surface.rotate_xyz(object, scene)
        render(scene, directory)

    elif option == 3:
        object, scene = surface.generate_obj_scene(vertex, faces)
        directory = "_spin_raw"
        surface.spin_bf(object, scene)
        render(scene, directory)

    elif option == 4:
        directory = "_multi_rotate_raw"
        object, object1, object2, scene = surface.generate_multi_obj_scene(
            vertex, faces)
        surface.multi_rotate(object, object1, object2, scene)
        render(scene, directory)

    elif option == 5:
        object, scene = surface.generate_obj_scene(vertex, faces)
        directory = "_rwalk_raw"
        surface.translate(object, scene)
        render(scene, directory)


# func_dict = {'smooth': smooth, 'raw': raw}
if __name__ == "__main__":
    choices = ["Raw", "Smooth"]
    choice = user_pick(choices)

    if choice == 1:
        raw()

    elif choice == 2:
        smooth()

    # command = input("Enter 'raw' or 'smooth': ")
    # func_dict[command]()
bpy.ops.wm.quit_blender()
bpy.ops.wm.window_close()
