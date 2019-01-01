""" test script to compare matplotlib plotter to glumpyplotter """
import bertini_real

#  bertini_real.plot.plot()

# color functions
def f1(x,y,z):
    return x**2 + y

def f2(x,y,z):
    return y**2 + z

def f3(x,y,z):
    return z**2 + x

function = f1, f2, f3

bertini_real.glumpyplotter.plot(function)

"""
these are just hardcoded examples of other color functions
"""

#  r = x**2 + y**2 + z**2
#  g = x + y**2 + z**2
#  b = x + y + z**2

# this is my fav
#  r = x**2
#  g = y**2
#  b = z**2
