import os
#def function_name(parameters)
#function doc string under
def parse_directory_name(directory_name = 'Dir_Name'):
	""" Parses file that contains the directory name, the MPtype, and the dimension
	 Returns a list [directory, MPtype, dimension]
	"""
	#os path manipulations. function returns true if path is an existing file
	if not os.path.isfile(directory_name):
		print "File does not exist, please run Bertini Real"
		return
	#open file, reads each line for the directory, mp type, and dimension and returns it in a list
	with open(directory_name, 'r') as f:
		directory = f.readline().replace('\n', '')
		MPtype = f.readline().replace('\n', '')
		dimension = f.readline().replace('\n', '')
	return [directory, MPtype, dimension]




#def defines function, name of function is parse_decomposition, parameter is the directory
def parse_decomposition(directory):
	"""
	Reads data from decomp file.
	Inputs: current directory

	Returns: List containing data to be stored in BRinfo class instance
		[Pi, Patch_Vectors, radius, center]

	"""
#checks if path file is a directory
#what does '/decomp' and '/r' do?
#reads for input file name, input variables and dimensions
#splits into two separate parts, separate integers
#turns them into integers
	if not os.path.isfile(directory + '/decomp'):
		print "did not find decomp at %s" %os.getcwd()
		return {}
	with open(directory + '/decomp', 'r') as f:
		inputFileName = f.readline().replace('\n','')
		num_variables_and_dimension = f.readline().replace('\n', '').split(' ')
		num_variables = int(num_variables_and_dimension[0])
		dimension = int(num_variables_and_dimension[1])

		Pi = [ [0,0] for i in xrange(num_variables - 1)]
#im kinda confused as to what this is doing
		for ii in xrange(dimension):
			numVars = f.readline()
			while numVars == '\n':
				numVars = f.readline()
			numVars = int(numVars.replace('\n', ''))
			for jj in xrange(numVars):
				Pi_Nums = f.readline().replace('\n', '').split(' ')
				if jj == 0:
					continue
				Pi[jj-1][ii] = complex(float(Pi_Nums[0]), float(Pi_Nums[1]))

#can we talk about what patches are?
#are patches extra blank lines?
		# Getting number of Patches and then patch data
		num_patches = f.readline()
		while num_patches == '\n':
			num_patches = f.readline()
		num_patches = int(num_patches.replace('\n', ''))


		Patch_Vectors = []
		for ii in xrange(num_patches):
			Patch_Vectors.append([])
			patch_size = f.readline()
			while patch_size == '\n':
				patch_size = f.readline()
			patch_size = int(patch_size.replace('\n', ''))
			for jj in xrange(patch_size):
				Patch_Vectors_data = f.readline().replace('\n','').split(' ')
				Patch_Vectors[ii].append(complex(float(Patch_Vectors_data[0]), float(Patch_Vectors_data[1])))
		#Get radius
		radius = f.readline()
		while radius == '\n':
			radius = f.readline()
		radius = radius.replace('\n', '').split(' ')
		radius = complex(float(radius[0]),float(radius[1]))

		centerSize = f.readline()
		while centerSize == '\n':
			centerSize = f.readline()
		centerSize = int(centerSize.replace('\n', ''))
		center = []
		for ii in xrange(centerSize):
			center_data = f.readline().replace('\n', '').split(' ')
			center.append(complex(float(center_data[0]), float(center_data[1])))
		return {'input file name':inputFileName,'Pi info': Pi, 'Patch Vectors':Patch_Vectors, "radius": radius, "center":center, "num patches":num_patches}

#reads data from surf file, returns number of faces, edges, midpoint slices, critical point singular_curve_multiplicites
#singular curves, and multiplicities
def parse_Surf(directory):
	""" Reads data from S.Surf file
	Inputs: current directory

	"""

	if not os.path.isfile(directory + '/S.surf'):
		print "S.surf does not exist in current directory: %s" %os.getcwd()
		return
	with open(directory + "/S.surf", 'r') as f:
		data = f.readline().replace('\n', '').split(' ')
		num_faces = int(data[0])
		num_edges = int(data[1])
		num_midpoint_slices = int(data[2])
		num_critical_point_slices = int(data[3])
		num_singular_curves = f.readline()
		while num_singular_curves == '\n':
			num_singular_curves = f.readline()
		num_singular_curves = int(num_singular_curves.replace('\n', ''))
		singular_curve_multiplicites = [ [0,0] for i in xrange(num_singular_curves)]
		multiplicites = f.readline().replace('\n', '').split(' ')
		index = 0
		for ii in xrange(num_singular_curves):
			for jj in xrange(2):
				singular_curve_multiplicites[ii][jj] = int(multiplicites[index])
				index += 1

	return [num_faces, num_edges, num_midpoint_slices, num_critical_point_slices,num_singular_curves, singular_curve_multiplicites]




def parse_Faces(directory):
	""" Reads Faces data from F.faces
	Inputs: current directory
	Returns: list with each element being a dictionary containing the face data
		Keys for each dictionary:
			"midpoint", "middle slice index", "top", "bottom"
			"system top", "system bottom", "num left". "left"
			"num right", "right"

	"""
	if not os.path.isfile(directory + '/F.faces'):
		print "F.faces file not found in current directory: %s" %os.getcwd()
		return
	with open(directory + '/F.faces') as f:
		num_faces = int(f.readline().replace('\n', ''))
		faces = [{} for i in xrange(num_faces)]
		for ii in xrange(num_faces):
			midpoint = f.readline()
			while midpoint == '\n':
				midpoint = f.readline()
			faces[ii]['midpoint'] = int(midpoint.replace('\n', ''))
			faces[ii]['middle slice index'] = int(f.readline().replace('\n', ''))

			top_bottom_data = f.readline().replace('\n', '').split(' ')
			faces[ii]['top'] = int(top_bottom_data[0])
			faces[ii]['bottom'] = int(top_bottom_data[1])

			system_data = f.readline().replace('\n', '').split(' ')
			faces[ii]['system top'] = system_data[0]
			faces[ii]['system bottom'] = system_data[1]

			faces[ii]['num left'] = int(f.readline().replace('\n', ''))
			if faces[ii]['num left'] == 0:
				faces[ii]['left'] = None
				f.readline()
			else:
				faces[ii]['left'] = [int(i) for i in f.readline().replace(' \n', '').split(' ')]

			faces[ii]['num right'] = int(f.readline().replace('\n', ''))
			if faces[ii]['num right'] == 0:
				faces[ii]['right'] = None
				f.readline()
			else:
				faces[ii]['right'] = [int(i) for i in f.readline().replace(' \n', '').split(' ')]

	return faces


def parse_Edges(directory):
	if not os.path.isfile(directory + '/E.edge'):
		print "E.edge file not found in current directory: %s" %os.getcwd()
		return {'number of edges': 0, 'edges': []}

	with open(directory + '/E.edge', 'r') as f:
		curves = {}
		curves['number of edges'] = int(f.readline().replace('\n', ''))
		curves['edges'] = [ [0,0,0] for i in xrange(curves['number of edges'])]
		for ii in xrange(curves['number of edges']):
			edges = f.readline()
			while edges == '\n':
				edges = f.readline()
			edges = edges.replace(' \n', '').split(' ')
			for jj in range(3):
				curves['edges'][ii][jj] = int(edges[jj]) + 1
				

	return curves



def parse_Curve_Sampler(directory):
	filename = directory + '/samp.curvesamp'
	if not os.path.isfile(filename):
		return [[]] #no curve sampling to parse, so return empty sampling

	with open(filename, 'r') as f:
		num_edges = int(f.readline().replace('\n', ''))
		f.readline() # read blank line.
		sampler_data = []

		for ii in xrange(num_edges):
			num_samples = int(f.readline().replace('\n', ''))
			temp = []
			thing = f.readline().replace('\n', '').split()
			for jj in thing:
				temp.append(int(jj))
			sampler_data.append(temp)
			f.readline() # read blank line.

		return sampler_data



def parse_Surface_Sampler(directory):
	filename = directory + '/samp.surfsamp'
	if not os.path.isfile(filename):
		return []
	else:
		pass
		# Finish this when running an example where a surface sample is produced
