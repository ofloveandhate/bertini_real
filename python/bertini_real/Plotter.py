
#Nicolle Ho
#University of Notre Dame
#Spring 2017

import os
from BRdata import BRdata
from Surface import Surface, Curve
import Util
import dill
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class BRplotter(object):

    def __init__(self, data = []):
        self.decomposition = data


    def plot(self):
        print("plotting object of dimension " + str(self.decomposition.dimension))
        self.fig = plt.figure()

        if self.decomposition.num_variables==2:
            self.ax = self.fig.add_subplot(111)
        else:
            self.ax = self.fig.add_subplot(111, projection='3d')
            self.ax.set_zlabel("z")
        self.PlotVertices()

        if self.decomposition.dimension == 1:
            self.PlotCurve(self.decomposition)
        elif self.decomposition.dimension == 2:
            self.PlotSurface(self.decomposition)

        plt.xlabel("x")
        plt.ylabel("y")
        #plt.xlabel
        plt.title(os.getcwd().split(os.sep)[-1])
        plt.show()




    def ReadMostRecent(self):
        filenum = Util.highest_filenumber()

        fileName = "BRdata" + str(filenum) + ".pkl"

        print("reading from file " + fileName)

        fileObject = open(fileName,'rb')
        self.decomposition = dill.load(fileObject)
        fileObject.close()

    def PlotVertices(self):

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
            self.ax.scatter(xs, ys, zs, zdir='z', s=20, c=None, depthshade=True)#v['point'][

    def PlotCurve(self, curve):
        num_nondegen_edges=0
        nondegen=[]
        for i in range(curve.curve.num_edges):
            e=curve.curve.edges[i]
            if e[0]!=e[1]!=e[2]:
                num_nondegen_edges=num_nondegen_edges+1
                nondegen.append(i)

        colormap = plt.cm.inferno
        color_list=([colormap(i) for i in np.linspace(0, 1,num_nondegen_edges)])

        for i in range(num_nondegen_edges):
            color=color_list[i]
            self.PlotEdge(curve,nondegen[i],color)

    def PlotEdge(self, curve, edge_index,color):

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


    def PlotSurface(self, surf):
        print("PlotSurface unimplemented yet")





def plot(data):
    b = BRplotter(data)
    b.plot();




if __name__ == "__main__":
     b = BRplotter()
     b.ReadMostRecent()
     b.plot()
