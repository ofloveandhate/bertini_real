#Nicolle Ho
#University of Notre Dame
#Spring 2017

# Danielle Brake
# Fall 2018

import os
from bertini_real.data import BRData
from bertini_real.surface import Surface, Curve
import bertini_real.util
import dill
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection



class StyleOptions(object):
	def __init__(self):
		self.line_thickness = 2 # there is no code using this yet.  write it.
		self.colormap = plt.cm.inferno

class VisibilityOptions(object):
	def __init__(self):
		self.vertices = True
		self.samples = True
		self.raw = True

		self.autorawsamples = True

		self.labels = False


class Options(object):
	def __init__(self):
		self.style = StyleOptions()
		self.visibility = VisibilityOptions()


class Plotter(object):
	def __init__(self, data = None, options = Options()):

		if data is None:
			self.decomposition = bertini_real.data.ReadMostRecent()
		else:
			self.decomposition = data

		self.options = options



	def plot(self):

		print("plotting object of dimension " + str(self.decomposition.dimension))
		
		self.make_figure()
		self.make_axes()
			
		self.main()
		

		self.label_axes()
		self.apply_title()


		plt.show()



	def main(self):

		self.points = self.extract_points()

		if self.options.visibility.vertices:
			self.plot_vertices()

		if self.decomposition.dimension == 1:
			self.plot_curve()
		elif self.decomposition.dimension == 2:
			self.plot_surface()


		

	def make_figure(self):
		self.fig = plt.figure()


	def make_axes(self):
		if self.decomposition.num_variables==2:
			self.ax = self.fig.add_subplot(1,1,1)
		else:
			self.ax = self.fig.add_subplot(1,1,1, projection='3d')


	def apply_title(self):
		plt.title(os.getcwd().split(os.sep)[-1])


	def label_axes(self):
		#todo: these should be set from the decomposition, not assumed to be x,y,z
		self.ax.set_xlabel("x")
		self.ax.set_ylabel("y")
		if self.decomposition.dimension is 2:
			self.ax.set_zlabel("z")

	





	'''
	renders all vertices

	todo: make them colored based on a function
	'''
	def plot_vertices(self):

		# refacotred version
		xs,ys,zs = self.make_xyz()


		if self.decomposition.num_variables==2:
			self.ax.scatter(xs, ys)
		else:
			self.ax.scatter(xs, ys, zs,
							zdir='z', s=5, c=None, depthshade=True)


	# this works well for plot_vertices
	# how can we make it work well for all the other methods???
	# returns 3 separate lists
	def make_xyz(self):
		xs = []
		ys = []
		zs = []

		for v in self.decomposition.vertices:
			xs.append(v['point'][0].real)
			ys.append(v['point'][1].real)
			if self.decomposition.num_variables>2:	
				zs.append(v['point'][2].real)

		return xs,ys,zs



	# would this be a better implementation than make_xyz???
	# returns one list
	def extract_points(self):
		points = []

		# get this from plot_samples.py
		for v in self.decomposition.vertices:
			#allocate 3 buckets to q
			q=[None]*3

			for i in range(3):
				#q[0],q[1],q[2]
				q[i]=v['point'][i].real
			points.append(q)

		return points




	def plot_curve(self):

		curve = self.decomposition.curve # a local unpacking

		should_plot_raw = self.options.visibility.raw
		should_plot_samp = self.options.visibility.samples and curve.sampler_data is not None
		if self.options.visibility.autorawsamples:
			if should_plot_samp:
				should_plot_raw = False;

		self.determine_nondegen_edges()

		if should_plot_raw:
			self.plot_raw_edges()

		
		if should_plot_samp:
			self.plot_edge_samples()
	

	def plot_raw_edges(self):
		curve = self.decomposition.curve # a local unpacking

		num_nondegen_edges = len(self.nondegen)

		colormap = self.options.style.colormap
		color_list=[colormap(i) for i in np.linspace(0, 1,num_nondegen_edges)]

		# instead of v['point']...etc look up into "self.points"
		for i in range(num_nondegen_edges):
			color=color_list[i]
			edge_index = self.nondegen[i]
			xs=[]
			ys=[]
			zs=[]
			inds = curve.edges[edge_index]
			for i in inds:
				v = self.decomposition.vertices[i]
				xs.append(v['point'][0].real)
				ys.append(v['point'][1].real)
				if self.decomposition.num_variables>2:
					zs.append(v['point'][2].real)

			if self.decomposition.num_variables==2:
				self.ax.plot(xs, ys, c=color)#v['point'][
			else:
				self.ax.plot(xs, ys, zs, zdir='z', c=color)#v['point']


	def plot_edge_samples(self):

		num_nondegen_edges = len(self.nondegen)

		colormap = self.options.style.colormap
		color_list=[colormap(i) for i in np.linspace(0, 1,num_nondegen_edges)]

		for i in range(num_nondegen_edges):
			color=color_list[i]
			edge_index = self.nondegen[i]
			xs=[]
			ys=[]
			zs=[]
			inds = self.decomposition.curve.sampler_data[edge_index]
			for i in inds:
				v = self.decomposition.vertices[i]
				xs.append(v['point'][0].real)
				ys.append(v['point'][1].real)
				if self.decomposition.num_variables>2:
					zs.append(v['point'][2].real)

			if self.decomposition.num_variables==2:
				self.ax.plot(xs, ys, c=color)#v['point'][
			else:
				self.ax.plot(xs, ys, zs, zdir='z', c=color)#v['point']
	
	def determine_nondegen_edges(self):
		curve = self.decomposition.curve # a local unpacking
		self.nondegen=[]
		for i in range(curve.num_edges):
			e=curve.edges[i]
			if e[0]!=e[1]!=e[2]:
				self.nondegen.append(i)

	



	def plot_surface(self):
		surf = self.decomposition # a local unpacking

		if self.options.visibility.samples:
			print("testing plot vertices")
			self.plot_surface_samples()

		if self.options.visibility.raw:
			self.plot_surface_raw()


	def plot_surface_samples(self):
		points = self.points

		tuples = self.decomposition.surface.surface_sampler_data

		colormap = self.options.style.colormap
		color_list=[colormap(i) for i in np.linspace(0, 1,num_nondegen_edges)]

		T = []

		for i in range(len(tuples)):

			# Initialize T here

			for tri in tuples[i]:
				f = int(tri[0])
				s = int(tri[1])
				t = int(tri[2])

				k = [points[f],points[s],points[t]]
				T.append(k)

		# add the collection here, with colors

		self.ax.add_collection3d(Poly3DCollection(T))



	def plot_surface_raw(self):
		print("plot_surface_raw unimplemented")


def plot(data = None, options = Options()):

	b = Plotter(data, options=options)
	b.plot();
	return b




