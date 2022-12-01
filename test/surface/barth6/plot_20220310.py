import bertini_real
import matplotlib.pyplot as plt

# bertini_real.data.gather_and_save()



decomposition = bertini_real.data.read_most_recent()

print(f'done reading data')

pieces = decomposition.separate_into_nonsingular_pieces()
print(f'separated into {len(pieces)} pieces')




if False:
	plotter = bertini_real.plot.Plotter()
	plotter.options.render.vertices = False
	plotter.options.render.surface_raw = False

	plotter.options.style.colormap = plt.cm.nipy_spectral

	# plotter.plot(decomposition)
	plotter.plot(pieces)



# plotter = bertini_real.plot.plot(pieces[0])


for p in pieces:
	p.solidify_smooth(distance=0.01)
	p.solidify_raw(distance=0.01)
