import os
from BRdata import BRdata
import Util

class BRplotter(object):

    def __init__(self):
        self.decomposition = []


    def plot(self):
        print("plotting object of dimension " + str(self.decomposition.dimension))


    def ReadMostRecent(self):
        import pickle

        filenum = Util.highest_filenumber()

        fileName = "BRdata" + str(filenum) + ".pkl"
        
        print("reading from file " + fileName)

        fileObject = open(fileName,'rb')
        self.decomposition = pickle.load(fileObject)
        fileObject.close()

        



if __name__ == "__main__":

     b = BRplotter()

     b.ReadMostRecent()

     b.plot()

