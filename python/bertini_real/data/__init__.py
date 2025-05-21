"""
    :platform: Unix, Windows
    :synopsis: The data module contains methods to find diretcory, gather vertices,read msot recent, gather and gather & save.

"""


from __future__ import print_function
import os
import bertini_real.parse as parse
from bertini_real.dehomogenize import dehomogenize
from bertini_real.surface import Surface, Curve
import bertini_real.util
from bertini_real.vertex import Vertex
import dill
import numpy as np



def find_directory(directory_path):
    """ Find the directory of decomposition data based on directory path

        :param directory_path: The path of directory

    """
    directories = directory_path.split('/')
    if len(directories) > 1:
        directory = directories[-1]
    else:
        directory = directory_path

    return directory

def gather_vertices(directory):
    """ Gather vertices

        :param directory: The path of directory

    """
    vertex_file_name = "V.vertex"
    if os.path.isfile("%s/V_samp.vertex" % directory):
        vertex_file_name = "V_samp.vertex"

    with open("%s/%s" % (directory, vertex_file_name), 'r') as f:
        # read first line and get number of vertices, number of projections,
        # number of natural vars, and number of file names
        line = f.readline().split(' ')
        num_vertices = int(line[0])
        num_projections = int(line[1])
        num_natural_vars_incl_hom_coord = int(line[2])
        num_variables = num_natural_vars_incl_hom_coord - 1
        num_filenames = int(line[3].replace('\n', ''))
        filenames = []

        # Skips data not used
        for i in range(0, (num_projections * num_natural_vars_incl_hom_coord)):
            skip_this_line = f.readline()
            while skip_this_line == '\n':
                skip_this_line = f.readline()
        _ = f.readline()
        if line == '\n':
            _ = f.readline()

        # Gets file names and stores them in filenames
        for ii in range(num_filenames):
            _ = f.readline()
            filenames.append(f.readline().replace('\n', ''))

        vertices = [None for i in range(num_vertices)]

        for ii in range(num_vertices):
            line = f.readline()

            while line == '\n':
                line = f.readline()
            number_of_variables = int(line)

            temporary_point = []
            for jj in range(number_of_variables):
                complex_num = f.readline().split(' ')
                real_part = float(complex_num[0])
                imaginary_part = float(complex_num[1])

                temporary_point.append(complex(real_part, imaginary_part))

            x = dehomogenize(temporary_point[0:num_natural_vars_incl_hom_coord])
            point = np.array(x)
            line = f.readline()

            while line == '\n':
                line = f.readline().replace('\n', '')

            num_projection_values = int(line)
            proj = []
            for ll in range(num_projection_values):
                complex_num = f.readline().split(' ')
                real_part = float(complex_num[0])
                imaginary_part = float(complex_num[1])
                proj.append(
                    complex(real_part,
                            imaginary_part))

            line = f.readline().replace('\n', '')

            while line == '\n':
                line = f.readline()

            filename_index = int(line.replace('\n', ''))


            all_filename_indices = "feature added in 1.9, your data is too old.  run a new decomposition using 1.9 or higher"
            if bertini_real.__version_info__ >= (1,9):
                line = f.readline()
                num_filename_indices = int(line)
                line = f.readline()
                all_filename_indices = [int(n) for n in line.split()]
                if len(all_filename_indices) != num_filename_indices:
                    raise RuntimeError(f"incorrect number of filename indices read, expected {num_filename_indices} but got {len(all_filename_indices)}")
            
            line = f.readline()

            vertextype = int(line.replace('\n', ''))

            path_numbers_ending_here = "feature added in 1.8, your data is too old.  run a new decomposition using 1.8 or higher"
            if bertini_real.__version_info__ >= (1,8):
                line = f.readline()
                num_paths_ending_here = int(line)
                line = f.readline()
                path_numbers_ending_here = [int(n) for n in line.split()]

            vertices[ii] = Vertex(point, filename_index, proj, vertextype, path_numbers_ending_here, all_filename_indices)

        return vertices, filenames

    def autosave(self):
        """ Going to remove this method """
        fn = util.next_filenumber()
        filename = "BRdata" + str(fn) + ".pkl"
        fileObject = open(filename, 'wb')
        dill.dump(b, fileObject)
        fileObject.close()
        return filename

def read(filename):
    """
    a nifty function for using dill to load a decomposition from a pickle file, given its name
    """

    print("reading from file " + filename)

    fileObject = open(filename, 'rb')
    decomposition = dill.load(fileObject)
    fileObject.close()

    return decomposition

def read_most_recent():
    """ Reads the most recent decomposition, and returns it."""
    filenum = bertini_real.util.highest_filenumber()

    filename = "BRdata" + str(filenum) + ".pkl"

    return read(filename)


def gather():
    """ 
    reads a decomposition from raw bertini_real output into memory, BUT NOT TO DISK.  use `gather_and_save` if you want to combine the steps.
    """
    
    directory_info = parse.parse_directory_name()
    directory = find_directory(directory_info[0])

    print("gathering data from " + directory)
    dimension = int(directory_info[2])
    if dimension == 1:
        # polynomial is a curve
        return Curve(directory, False)
    elif dimension == 2:
        # polynomial is a surface
        return Surface(directory, False)

def gather_and_save():
    """ Gather and save data """
    a = bertini_real.util.next_filenumber()

    filename = "BRdata" + str(a) + ".pkl"
    fileObject = open(filename, 'wb')
    b = gather()

    print("saving to file " + filename)

    import dill
    dill.dump(b, fileObject)
    fileObject.close()
