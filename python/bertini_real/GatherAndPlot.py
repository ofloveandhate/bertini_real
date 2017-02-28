import bertini_real as br

from bertini_real.BRdata import BRdata
from bertini_real.Surface import Surface, Curve
import bertini_real.Plotter
import bertini_real.Util as Util

a = Util.next_filenumber()

fileName = "BRdata" + str(a) + ".pkl"
fileObject = open(fileName,'wb')

b = BRdata()

print("saving to file " + fileName)

import dill
dill.dump(b,fileObject)
fileObject.close()



p = br.Plotter.BRplotter()

p.ReadMostRecent()

p.plot()
