import os
import bertini_real


def rotate_xyz():
    script_path = os.path.dirname(
        bertini_real.__file__) + "/anaglypy/rotate_xyz.py"
    os.system("blender -P " + script_path)


def rotate_z():
    script_path = os.path.dirname(
        bertini_real.__file__) + "/anaglypy/rotate_z.py"
    os.system("blender -P " + script_path)


if __name__ == "__main__":
    rotate_xyz()
    rotate_z()
