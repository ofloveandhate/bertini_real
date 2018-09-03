
import bertini_real.brdata
import bertini_real.surface
import bertini_real.curve
import bertini_real.util




def gather():

    a = util.next_filenumber()

    fileName = "BRdata" + str(a) + ".pkl"
    fileObject = open(fileName,'wb')
#file parsing
    b = brdata.BRData()

    print("saving to file " + fileName)

    import dill
    #b is written to fileObject through dump function
    dill.dump(b,fileObject)
    #closes the file
    fileObject.close()

def gather_and_plot():
	gather()
	bertini_real.Plot.plot()