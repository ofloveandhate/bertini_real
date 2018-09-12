import os
from bertini_real.data import BRData
from bertini_real.surface import Surface, Curve
import bertini_real.util
import dill
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri

import bertini_real

fig = plt.figure(figsize=plt.figaspect(0.5))

sphere_data = bertini_real.data.ReadMostRecent()
sphere_tuples = sphere_data.surface.surface_sampler_data

x = []
y = []
z = []

for i in sphere_tuples:
    x.append(i[0])
    y.append(i[1])
    z.append(i[2])

# tri = mtri.Triangulation(x,y)


ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.plot_trisurf(x,y,z)

plt.show()
