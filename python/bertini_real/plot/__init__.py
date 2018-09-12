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
			self.ax = self.fig.add_subplot(111)
		else:
			self.ax = self.fig.add_subplot(111, projection='3d')


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

		xs = []
		ys = []
		zs = []

		for v in self.decomposition.vertices:
			xs.append(v['point'][0].real)
			ys.append(v['point'][1].real)

			if self.decomposition.num_variables>2:
				zs.append(v['point'][2].real)


		if self.decomposition.num_variables==2:
			self.ax.scatter(xs, ys)#v['point'][
		else:
			self.ax.scatter(xs, ys, zs, zdir='z', s=5, c=None, depthshade=True)#v['point'][






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

	



	# what exactly is this supposed to do
	def plot_surface(self):
		surf = self.decomposition # a local unpacking
		print("plot_surface unimplemented yet.")

		if self.options.visibility.samples:
			self.PlotSurfaceSamples()

		if self.options.visibility.raw:
			self.plot_surface_raw()


	def plot_surface_samples(self):
		print("plot_surface_samples unimplemented yet.")


	def plot_surface_raw(self):
		print("plot_surface_raw unimplemented")


def plot(data = None, options = Options()):

	b = Plotter(data, options=options)
	b.plot();
	return b




