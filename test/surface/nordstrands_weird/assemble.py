import bertini_real as br
import matplotlib.pyplot as plt
import dill
# br.data.gather_and_save()

surf = br.data.read_most_recent()

# fileObject = open('BRdata3.pkl', 'rb')
# surf = dill.load(fileObject)
# fileObject.close()

surf.write_piece_data()

pieces = surf.separate_into_nonsingular_pieces()


br.surface.copy_all_scad_files_here()
# 
# plotter = br.plot.plot(pieces)

# plotter = br.plot.plot(pieces[0])

# surf.export_raw()
# surf.export_smooth()

for p in pieces:
	# p.export_smooth()
	p.solidify_smooth(0.02)

	print(p.point_singularities() )
# 	p.solidify_raw(0.02)



plotter = br.plot.Plotter()
plotter.options.render.vertices = False
plotter.options.render.surface_raw = False

plotter.options.style.colormap = plt.cm.nipy_spectral

# plotter.plot(decomposition)
plotter.plot(pieces)

