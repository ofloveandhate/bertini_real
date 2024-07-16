

import bertini_real as br

br.data.gather() # do this once after decomposing and sampling the surface. If 

surface = br.data.read_most_recent()

pieces = surface.separate_into_nonsingular_pieces()

surface.write_piece_data()

for p in pieces:
  p.export_smooth()
