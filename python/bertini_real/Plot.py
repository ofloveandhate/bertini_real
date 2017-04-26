import bertini_real.Plotter
import bertini_real as br

def plot():

    p = br.Plotter.BRplotter()
    #p takes the action ReadMostRecent()
    p.ReadMostRecent()
    #member functions have to be called as a member function
    p.plot()


if __name__ == "__main__":
    plot()
