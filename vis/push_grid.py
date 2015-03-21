import sys
import json
import math
import numpy as np
from matplotlib.collections import Collection
from matplotlib.artist import allow_rasterization
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import matplotlib as mpl

#mpl.rcParams['axes.linewidth'] = .75
mpl.rcParams['text.usetex'] = True

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

def insert_rasterized_contour_plot(c, ax = None):
    collections = c.collections
    for _c in collections:
        _c.remove()
    cc = ListCollection(collections, rasterized=True)
    if ax is None :
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

mpl.rcParams['figure.figsize'] = 6, 3

#fig, (ax1, ax2) = plt.subplots(1, 2)


l_y = []
m_x = []
c_z = []

for i in range(shape[0]) :
    for j in range(shape[1]) :
        l_y.append(-Y[i, j]*math.sin(math.pi*X[i, j]/180.0))
        m_x.append(Y[i, j]*math.cos(math.pi*X[i, j]/180))


l_y = np.reshape(l_y, (shape[0], shape[1]))
m_x = np.reshape(m_x, (shape[0], shape[1]))

fig, (ax1, ax2) = plt.subplots(1, 2)
cs1 = ax1.contourf(X, Y, Z, levels=levels, cmap='Greys')
insert_rasterized_contour_plot(cs1, ax1)

ax1.invert_xaxis()
plt.setp(ax1, xticks=[45, 90, 135], yticks=[2, 6, 10, 14, 18,  22])
ax1.tick_params(axis='both', which='major', labelsize=10)
ax1.set_xlabel('$\\theta$', fontsize=12)
ax1.set_ylabel('$\sqrt{\lambda^2 + \mu^2}$', fontsize=12)
ax1.axis('tight')

yticks = ax1.yaxis.get_major_ticks()
yticks[0].tick1On = False
yticks[0].tick2On = False
yticks[-1].tick1On = False
yticks[-1].tick2On = False

cs2 = ax2.contourf(m_x, l_y, Z, levels=levels, cmap='Greys')
insert_rasterized_contour_plot(cs2, ax2)

ax2.tick_params(axis='both', which='major', labelsize=10)
ax2.set_xlabel('$\mu$', fontsize=12)
ax2.set_ylabel('$\lambda$', fontsize=12)
ax2.axis('tight')

plt.tight_layout()

cbar = fig.colorbar(cs2, ax = [ax1, ax2], ticks=[0, 16, 32, 48, 64, 80, 96], fraction=.1)
cbar.solids.set_edgecolor('face')

cbar.ax.tick_params(labelsize=10)
ycticks = cbar.ax.yaxis.get_major_ticks()
ycticks[0].tick2On = False
ycticks[-1].tick2On = False

plt.draw()

plt.savefig('temp_push.eps', format='eps', dpi=300)
plt.show()
