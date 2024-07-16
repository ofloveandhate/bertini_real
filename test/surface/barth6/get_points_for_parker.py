import bertini_real
import matplotlib.pyplot as plt

import numpy as np









# only do this next line if need to gather.  it makes a new pkl file.
# bertini_real.data.gather_and_save()



decomposition = bertini_real.data.read_most_recent()

print(f'done reading data')

pieces = decomposition.separate_into_nonsingular_pieces()
print(f'separated into {len(pieces)} pieces')

num_variables = decomposition.num_variables

# print(decomposition.sampler_data)

all_singularities_indices = []

a_compact_piece = None

for p in pieces:
	if p.is_compact() and not a_compact_piece:
		a_compact_piece = p
	print(p.point_singularities())
	all_singularities_indices.extend(p.point_singularities())


s = [decomposition.vertices[ind].point[:num_variables].real for ind in set(all_singularities_indices)]

all_singularities = np.array([p for p in s if np.sqrt(np.sum(p**2))<3.1])
print(all_singularities.shape)


all_singularities.tofile('barth6_singularities.csv', sep = ',')








all_vertices = np.array( [v.point[:num_variables].real for v in decomposition.vertices] )
all_vertices.tofile('all_points.csv', sep = ',')



fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(all_singularities[:,0],all_singularities[:,1],all_singularities[:,2])
plt.show()


# points_on_one_compact_cone = p.face_points()
# points_on_one_compact_cone.tofile('barth6_one_cone.json',sep=',')







