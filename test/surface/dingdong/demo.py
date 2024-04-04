import bertini_real as br

surface = br.data.read_most_recent()

pieces = surface.separate_into_nonsingular_pieces()

surface.write_piece_data()

br.surface.copy_all_scad_files_here()

for p in pieces:
	p.export_smooth()


br.plot.plot(pieces)