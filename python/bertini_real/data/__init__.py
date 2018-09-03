from __future__ import print_function
import os
import bertini_real.parse as parse

from bertini_real.dehomogenize import dehomogenize
from bertini_real.surface import Surface, Curve
import bertini_real.util
import dill
import numpy as np


'''
implemented in a way such that every decomposition is is a BRData,
and then if it's a curve, it gets a curve, and if its a surface, it gets a surface.

i'd rather it was the data type itself, and modeled is-a rather than has-a, but this is what there is right now.  

todo: re-implement so that is-a instead of has-a is the thing to do. 
all related code could be much more generic.

'''
class BRData(object):
    def __init__(self, autoload = True):
        self.filenames = []
        self.num_vertices = 0
        self.vertices = []
        if autoload:
            self.gather()

    def __str__(self):
        result = ""

        result += "decomposition of dimension {}".format(self.dimension)

        return result

    def __repr__(self):
        return str(self)

    def gather(self):
        self.directory_info = bertini_real.parse.parse_directory_name()
        self.find_directory(self.directory_info[0])

        print("gathering data from " + self.directory)
        self.dimension = int(self.directory_info[2])
        # gather vertices
        self.gather_vertices()
        if self.dimension == 1:
            # polynomial is a curve
            self.gather_curve(self.directory)
        elif self.dimension == 2:
            # polynomial is a surface
            self.gather_surface(self.directory)
        print("done gathering decomposition")



    def find_directory(self,directory_path):
        directories = directory_path.split('/')
        if len(directories) > 1:
            self.directory = directories[-1]
        else:
            self.directory = directory_path



    def gather_vertices(self):
        vertex_file_name = "V.vertex"
        # print self.directory
        if os.path.isfile("%s/V_samp.vertex" %self.directory):
            vertex_file_name = "V_samp.vertex"

        # open vertex file and read data from file in read mode
        with open("%s/%s" %(self.directory, vertex_file_name), 'r') as f:
            # read first line and get number of vertices, number of projections,
            # number of natural vars,and number of file names
            line = f.readline().split(' ')
            self.num_vertices = int(line[0])
            num_projections = int(line[1])
            num_natural_vars = int(line[2])
            self.num_variables = num_natural_vars-1;
            num_filenames = int(line[3].replace('\n', ''))

            # Skips data not used
            for i in range(0, (num_projections * num_natural_vars)):
                skip_this_line = f.readline()
                while skip_this_line == '\n':
                        skip_this_line = f.readline()
            skip_line = f.readline()
            if line == '\n':
                skip_line = f.readline()


            # Gets file names and stores them in self.filenames
            for ii in range(num_filenames):
                skip_line = f.readline()
                self.filenames.append(f.readline().replace('\n', ''))

            self.vertices = [ {} for i in range(self.num_vertices)]

            for ii in range(self.num_vertices):
                line = f.readline()

                while line == '\n':
                    line = f.readline()
                number_of_variables = int(line)
                self.vertices[ii]['point'] = []
                temporary_point = []
                for jj in range(number_of_variables):
                    complex_num = f.readline().split(' ')
                    real_part = float(complex_num[0])
                    imaginary_part = float(complex_num[1])

                    temporary_point.append(complex(real_part, imaginary_part))

                x = dehomogenize(temporary_point[0:num_natural_vars])
                self.vertices[ii]['point'] = np.array(x)
                line = f.readline()

                while line == '\n':
                    line = f.readline().replace('\n', '')


                num_projection_values = int(line)
                self.vertices[ii]['projection_value'] = []
                for ll in range(num_projection_values):
                    complex_num = f.readline().split(' ')
                    real_part = float(complex_num[0])
                    imaginary_part = float(complex_num[1])
                    self.vertices[ii]['projection_value'].append(complex(real_part, imaginary_part))

                line = f.readline().replace('\n', '')

                while line == '\n':
                    line = f.readline()

                self.vertices[ii]['input_filename_index'] = float(line.replace('\n', ''))
                line = f.readline()

                while line == '\n':
                    line = f.readline()

                self.vertices[ii]['type'] = float(line.replace('\n',''))
        return


    def gather_surface(self, directory):
        # creates new surface
        self.surface = Surface(directory)


    def gather_curve(self, directory):
        self.curve =  Curve(directory)

    def autosave(self):
        fn = util.next_filenumber()
        fileName = "BRdata" + str(fn) + ".pkl"
        fileObject = open(fileName,'wb')
        dill.dump(b,fileObject)
        fileObject.close()
        return fileName


'''
Reads the most recent decomposition, and returns it.
'''
def ReadMostRecent():
        filenum = bertini_real.util.highest_filenumber()

        fileName = "BRdata" + str(filenum) + ".pkl"

        print("reading from file " + fileName)

        fileObject = open(fileName,'rb')
        decomposition = dill.load(fileObject)
        fileObject.close()

        return decomposition



def gather():

    a = bertini_real.util.next_filenumber()

    fileName = "BRdata" + str(a) + ".pkl"
    fileObject = open(fileName,'wb')
#file parsing
    b = BRData()

    print("saving to file " + fileName)

    import dill
    #b is written to fileObject through dump function
    dill.dump(b,fileObject)
    #closes the file
    fileObject.close()

