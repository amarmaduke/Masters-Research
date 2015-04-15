import sys
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from  matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = .75
mpl.rcParams['text.usetex'] = True

lw = .4

## Script Inputs
#view = [-25, 150, -1, 50]
view = [-25, 150, -1, 20]

path = '/home/owner/dev/Masters-Research/sims/many/runs/reference/'
data = {}
for dir_entry in os.listdir(path) :
    dir_entry_path = os.path.join(path, dir_entry)
    if os.path.isfile(dir_entry_path) :
        with open(dir_entry_path, 'r') as my_file :
            data[dir_entry] = my_file.read()

f, axarray = plt.subplots(3, 3)
c = 1

f.set_size_inches(7,4,forward=True) # landscape
# f.set_size_inches(6,5,forward=True)

for d in data :
    json_data = json.loads(data[d])
    n = json_data['n']
    m = json_data['m']
    delta = json_data['delta']
    N = n*m

    sindex, t = 0, 0

    try :
        while True:
            index = 'sindex%dtq%d' % (sindex,t,)
            test = json_data[index]
            t = t + 1
    except KeyError:
        t = t - 1

    print('Processing... %d' % (c))
    index = 'sindex%dtq%d' % (sindex,t,)
    stuff = json_data[index]
    a = np.asarray(stuff)
    lamb = json_data['lambda']
    mu = json_data['mu']

    i = int((lamb - 5)/10)
    j = int(mu/5) + 1
    ax = plt.subplot(3, 3, 3*i + j)

    #ax.set_aspect('equal',adjustable='box')
    ax.set_xlim(view[0], view[1])
    ax.set_ylim(view[2], view[3])

    if mu == 0 :
        plt.ylabel('$\lambda = %d$' % lamb, fontsize=12)
    if lamb == 25 :
        plt.xlabel('$\mu = %d$' % mu, fontsize=12)

    plt.tick_params(axis='both', which='both', bottom='off', top='off',
        right='off', left='off', labelbottom='off', labelleft='off')

    for j in range(m) :
        x2, y2 = delta[j], 0
        for i in range(n) :
            x1, y1 = a[i + j*n], a[i + j*n + N]
            ax.plot([x1, x2], [y1, y2], color='k', linewidth=lw)
            x2, y2 = x1, y1

    top_substrate_step = json_data['sub_h']
    top_substrate_count = json_data['sub_count']
    bottom_substrate_step = json_data['osub_h']
    bottom_substrate_count = json_data['osub_count']
    bottom_substrate_x = json_data['osub']

    x2, y2 = a[2*N], a[2*N + 1]
    for j in range(top_substrate_count) :
        x1, y1 = x2 + top_substrate_step, y2
        ax.plot([x1, x2], [y1, y2], color='k', linewidth=lw)
        x2, y2 = x1, y1

    x2, y2 = bottom_substrate_x, 0
    for j in range(bottom_substrate_count) :
        x1, y1 = x2 + bottom_substrate_step, y2
        ax.plot([x1, x2], [y1, y2], color='k', linewidth=lw)
        x2, y2 = x1, y1

    c = c + 1

plt.tight_layout()
plt.draw()
plt.savefig('temp.eps', format='eps', dpi=300)
plt.show()

