import sys
import json
import numpy as np
from matplotlib.collections import Collection
from matplotlib.artist import allow_rasterization
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms

class ListCollection(Collection):
     def __init__(self, collections, **kwargs):
         Collection.__init__(self, **kwargs)
         self.set_collections(collections)
     def set_collections(self, collections):
         self._collections = collections
     def get_collections(self):
         return self._collections
     @allow_rasterization
     def draw(self, renderer):
         for _c in self._collections:
             _c.draw(renderer)

def insert_rasterized_contour_plot(c):
    collections = c.collections
    for _c in collections:
        _c.remove()
    cc = ListCollection(collections, rasterized=True)
    ax = plt.gca()
    ax.add_artist(cc)
    return cc

json_raw = sys.stdin.read()
data = json.loads(json_raw)
n = data['n']
m = data['m']
print(m, n)
delta = data['delta']
N = n*m

x_range = data['range']
y_range = data['range2']

shape = [int((x_range[2] - x_range[0])/x_range[1]),
         int((y_range[2] - y_range[0])/y_range[1])]

grid = np.asarray(data['grid'])

X = grid[:,0]
Y = np.sqrt(grid[:,1]**2 + grid[:,2]**2)
Z = grid[:,6]

X = np.reshape(X, (shape[0], shape[1]))
Y = np.reshape(Y, (shape[0], shape[1]))
Z = np.reshape(Z, (shape[0], shape[1]))

levels = range(0, 97, 1)

plt.figure(1, (5.5, 4.25))
ax = plt.gca()
# Anti-aliasing rendering bugs, white-lines and artifacts etc
# In order to avoid that we put redundant information in the
# figure to "help" the renderer get pixels right
# Your mileage may vary
cs = plt.contourf(X, Y, Z, levels=levels, cmap="Greys")
#plt.contourf(X, Y, Z, levels=levels, cmap="Greys")
#plt.contourf(X, Y, Z, levels=levels, cmap="Greys")

insert_rasterized_contour_plot(cs)

plt.xticks([45, 90, 135], fontsize=10)
plt.xlabel('$\\theta$', fontsize=16)
plt.yticks([2, 6, 10, 14, 18,  22], fontsize=10)
plt.ylabel('$\sqrt{\lambda^2 + \mu^2}$', fontsize=16)
plt.axis('tight')
cbar = plt.colorbar(ticks=[0, 16, 32, 48, 64, 80, 96], fraction=.1)
cbar.solids.set_edgecolor('face')

yticks = ax.yaxis.get_major_ticks()
yticks[0].tick1On = False
yticks[0].tick2On = False
yticks[-1].tick1On = False
yticks[-1].tick2On = False

cbar.ax.tick_params(labelsize=10)
ycticks = cbar.ax.yaxis.get_major_ticks()
ycticks[0].tick2On = False
ycticks[-1].tick2On = False

plt.draw()

ax = plt.gca()
ax.invert_xaxis()

plt.savefig('temp_push.eps', format='eps', dpi=300)
plt.show()