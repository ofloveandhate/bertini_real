import os
from BRdata import BRdata
from Surface import Surface, Curve
import Util
import dill

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class BRplotter(object):

    def __init__(self, data = []):
        self.decomposition = data


    def plot(self):
        print("plotting object of dimension " + str(self.decomposition.dimension))
        print("NH selfdecomp: " + str(self.decomposition))
        self.fig = plt.figure()

        if self.decomposition.num_variables==2:
            self.ax = self.fig.add_subplot(111)
        else:
            self.ax = self.fig.add_subplot(111, projection='3d')
        self.PlotVertices()


        if self.decomposition.dimension == 1:
            self.PlotCurve(self.decomposition)
        elif self.decomposition.dimension == 2:
            self.PlotSurface(self.decomposition)

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
        print "number of variables:", self.decomposition.num_variables

        for v in self.decomposition.vertices:
            xs.append(v['point'][0].real)
            ys.append(v['point'][1].real)

            if self.decomposition.num_variables>2:
                zs.append(v['point'][2].real)


        if self.decomposition.num_variables==2:
            self.ax.scatter(xs, ys)#v['point'][
        else:
            self.ax.scatter(xs, ys, zs, zdir='z', s=20, c=None, depthshade=True)#v['point'][
        #print xs, ys, zs
        #for item in
    def PlotCurve(self, curve):
        for i in range(curve.curve.num_edges):
            self.PlotEdge(curve, i)

    def PlotEdge(self, curve, edge_index):

        xs = []
        ys = []
        zs = []

        inds = self.decomposition.curve.sampler_data[edge_index]
        for i in inds:
            print self.decomposition.vertices
            v = self.decomposition.vertices[i]
            print v
            xs.append(v['point'][0].real)
            ys.append(v['point'][1].real)
            if self.decomposition.num_variables>2:
                zs.append(v['point'][2].real)

        if self.decomposition.num_variables==2:
            self.ax.plot(xs, ys, c=None)#v['point'][
        else:
            self.ax.plot(xs, ys, zs, zdir='z', c=None)#v['point'][



        # get indicies from self.decomposition.curve.curve_sampler_data[]

        #those indicies index into self.decomposition.vertices

    def PlotSurface(self, surf):
        print("PlotSurface unimplemented yet")





def plot(data):
    b = BRplotter(data)
    b.plot();




if __name__ == "__main__":
     b = BRplotter()
     b.ReadMostRecent()
     b.plot()
