from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
import numpy as np
import bertini_real

#create a new figure
fig = plt.figure ()

#describe the position of the su plot
ax = fig.add_subplot (1, 1, 1, projection = '3d', aspect = 1)

#read the most recent data?
data = bertini_real.data.ReadMostRecent();

#tuples -  store surface sampler data
tuples = data.surface.surface_sampler_data


"""Extract points from vertices"""
def extractPoints(data):
	points = []

	for v in data.vertices:
		#allocate 3 buckets to q
		q=[None]*3

		for i in range(3):
			#q[0],q[1],q[2]
			q[i]=v['point'][i].real
		points.append(q)
	return points

#points - store extracted points
points = extractPoints(data)

#create an empty array T
T = []
#T=[[points[f],points[s],points[t]]]

"""Questions"""
# len(tuples) - 'int' object is not iterable?
# how to get size of tuples and size of list in tuples?
for i in range(2):
	
	for j in range(1000):
		f = int(tuples[i][j][0])
		s = int(tuples[i][j][1])
		t = int(tuples[i][j][2])
		#print(f,s,t)
		k = [points[f],points[s],points[t]]
		T.append(k)


ax.add_collection3d(Poly3DCollection(T))

plt.show()

