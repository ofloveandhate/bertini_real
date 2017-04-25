import bertini_real as br

from bertini_real.BRdata import BRdata
from bertini_real.Surface import Surface, Curve
import bertini_real.Plotter
import bertini_real.Util as Util

def gather():

    a = Util.next_filenumber()

    fileName = "BRdata" + str(a) + ".pkl"
    fileObject = open(fileName,'wb')
#file parsing
    b = BRdata()

    print("saving to file " + fileName)

    import dill
    #b is written to fileObject through dump function
    dill.dump(b,fileObject)
    #closes the file
    fileObject.close()

if __name__ == "__main__":
    gather()
