from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import numpy as np
import bertini_real
fig = plt.figure ()
ax = fig.add_subplot (1, 1, 1, projection = '3d', aspect = 1)

# extract and flatten points from the vertices

# decomposition.vertices[ii]['point']
sphere_data = bertini_real.data.ReadMostRecent();
sphere_tuples = sphere_data.surface.surface_sampler_data

#once!!!  store the result in self (the plotter object)

#vertices = []
#x = []
#y = []
#z = []

#for jj in sphere_tuples
#	vertices.append(

#vertices.append([ 0.17770898,  0.72315927,  0.66742804])
#vertices.append([-0.65327074, -0.4196453 ,  0.63018661])
#vertices.append([ 0.65382635,  0.42081934, -0.62882604])
#vertices.append([-0.17907021, -0.72084723, -0.66956189])
#vertices.append([-0.73452809,  0.5495376 , -0.39809158])
#vertices.append([ 0.73451554, -0.55094017,  0.39617148])


#actual_data = [[],[],[]]
#happens in the plot_samples function

#actual_data[0][0] = x.vertices[0]['point'][0]
#actual_data[0][1] = x.vertices[0]['point'][1]
#actual_data[0][2]



#for ii in range(num_faces):
	#already have this data
	# sampler_data[face_index] = [[4,0,1],
	# 		[4,1,3],
	# 		[4,3,2],
	# 		[4,2,0],
	# 		[5,0,1],
	# 		[5,1,3],
	# 		[5,3,2],
	# 		[5,2,0],
	# ]

#	actual_data = []
#	for t in sampler_data[ii]:
#		actual_data.append([vertices[t[0]],vertices[t[1]],vertices[t[2]]])

f = int(sphere_tuples[0][0])
s = int(sphere_tuples[0][1])
t = int(sphere_tuples[0][2])

f1 = sphere_data.vertices[f]
s1 = sphere_data.vertices[s]
t1 = sphere_data.vertices[t]

#print(f1, s2)

fx= [f1['point'][0].real]
fy= [f1['point'][1].real]
fz= [f1['point'][2].real]

sx=[s1['point'][0].real]
sy=[s1['point'][1].real]
sz=[s1['point'][2].real]

tx = [t1['point'][0].real]
ty = [t1['point'][1].real]
tz = [t1['point'][2].real]
#print(fx)


#ax.add_collection3d (Poly3DCollection (actual_data))
ax.plot([[fx,fy,fz],[sx,sy,sz],[tx,ty,tz]], marker=11)

# at the very end of the plot call
plt.show ()


