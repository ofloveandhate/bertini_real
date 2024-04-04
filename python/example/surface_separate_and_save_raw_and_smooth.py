import bertini_real as br
import matplotlib.pyplot as plt

import pathlib

pkl_files = list(pathlib.Path('.').glob('*.pkl'))

if not pkl_files:
	print('found no pickled data, so trying to gather_and_save')
	br.data.gather_and_save()



thickness = 0.01 # this is in object dimensions.  

decomposition = br.data.read_most_recent()

print(f'done reading data')

pieces = decomposition.separate_into_nonsingular_pieces()
print(f'separated into {len(pieces)} pieces')





for p in pieces:
	p.solidify_smooth(distance=thickness)
	p.solidify_raw(distance=thickness)
