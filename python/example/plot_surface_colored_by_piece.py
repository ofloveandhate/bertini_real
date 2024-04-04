# an example of splitting an algebraic surface into pieces separated by singularities, 
# and plotting each with a solid color

import bertini_real as br
import matplotlib.pyplot as plt

import pathlib

pkl_files = list(pathlib.Path('.').glob('*.pkl'))

if not pkl_files:
	print('found no pickled data, so trying to gather_and_save')
	br.data.gather_and_save()


decomposition = br.data.read_most_recent()
print(f'done reading data')


pieces = decomposition.separate_into_nonsingular_pieces()
print(f'separated into {len(pieces)} pieces')


# make a plotter object
plotter = br.plot.Plotter()

# set some options
plotter.options.render.vertices = False
plotter.options.render.surface_raw = False

plotter.options.style.colormap = plt.cm.viridis

# plot it!
plotter.plot(pieces)


