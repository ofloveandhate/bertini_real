import bertini_real as br
import matplotlib.pyplot as plt
import dill

def plot_pieces(pieces):
	plotter = br.plot.Plotter()
	plotter.options.render.vertices = False
	plotter.options.render.surface_raw = False

	plotter.options.style.colormap = plt.cm.nipy_spectral

	plotter.plot(pieces)


br.data.gather_and_save()

surf = br.data.read_most_recent()
surf.write_piece_data()
pieces = surf.separate_into_nonsingular_pieces()
br.surface.copy_all_scad_files_here()

for p in pieces:
	p.export_smooth()

plot_pieces(pieces)

