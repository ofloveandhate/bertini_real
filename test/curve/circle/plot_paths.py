import bertini_real as br
import numpy as np


C = br.data.gather()




import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(projection='3d')





colors = plt.cm.viridis(np.linspace(0, 1, len(C.vertices)))



for ii in range(len(C.vertices)):

	v = C.vertices[ii]

	path_nums = v.path_numbers_ending_here
	ax.scatter(np.real(v.point[0]), np.real(v.point[1]), 0, color=colors[ii,:])

	for n in np.unique(path_nums):
		p = br.paths.Path(f"paths/path_{n+1}")

		s = p.space[:,0:C.num_variables]
		ax.plot(np.real(s[:,0]), np.real(s[:,1]), np.abs(np.imag(s[:,1]))+np.abs(np.imag(s[:,0])),color=colors[ii,:])


plt.show()