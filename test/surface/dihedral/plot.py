import bertini_real as  br

dihedral = br.data.read('BRdata0.pkl')


options = br.plot.Options()

options.style.colormode = br.plot.ColorMode.BY_FUNCTION
options.style.color_function = lambda x: [x[0],x[1],x[2],0.75]
options.style.autotitle = False

options.render.curve_samples = True
options.render.defer_show = True

options.render.surface_curves=True

plotter = br.plot.Plotter(options=options)


plotter.plot(dihedral)

ax = plotter.ax

ax.set_xlim([-1.7, 1.7])
ax.set_ylim([-1.7, 1.7])
ax.set_zlim([-1.1, 1.1])

ax.axis('off')
import matplotlib.pyplot as plt
plt.tight_layout()

plotter.show()

