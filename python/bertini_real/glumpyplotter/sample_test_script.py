import bertini_real

def function(x,y,z):
    return x + y + z

fn = function

bertini_real.glumpyplotter.plot()
#  bertini_real.glumpyplotter.plot(cmap='jet')
#  bertini_real.glumpyplotter.plot(cmap='inferno', color_function=fn)
#  bertini_real.glumpyplotter.plot(cmap='inferno', color_function=fn, critical_curve=True)
