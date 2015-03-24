import sys
import json
import numpy as np
from matplotlib.collections import Collection
from matplotlib.artist import allow_rasterization
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = .75
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

grid = np.asarray(data['grid'])

X = grid[:,0]
Y = grid[:,1]
Z = grid[:,5]

X = np.reshape(X, (100, 100))
Y = np.reshape(Y, (100, 100))
Z = np.reshape(Z, (100, 100))

levels = range(0, 97, 1)

plt.figure(1, (3.5, 3))
ax = plt.gca()
cs = plt.contourf(X, Y, Z, levels=levels, cmap="Greys")

insert_rasterized_contour_plot(cs)

plt.xticks([0, 2, 4, 6, 8, 10], fontsize=10)
plt.xlabel('$\\beta$', fontsize=12)
plt.yticks([2, 4, 6, 8, 10], fontsize=10)
plt.ylabel('$\\varepsilon_-$', fontsize=12)
plt.axis('tight')

plt.tight_layout()

cbar = plt.colorbar(ticks=[0, 16, 32, 48, 64, 80, 96], fraction=.1)
cbar.solids.set_edgecolor('face')

yticks = ax.yaxis.get_major_ticks()
yticks[-1].tick1On = False
yticks[-1].tick2On = False

xticks = ax.xaxis.get_major_ticks()
xticks[-1].tick1On = False
xticks[-1].tick2On = False

cbar.ax.tick_params(labelsize=10)
ycticks = cbar.ax.yaxis.get_major_ticks()
ycticks[0].tick2On = False
ycticks[-1].tick2On = False

plt.savefig('temp_fs.eps', format='eps', dpi=300)
plt.show()
