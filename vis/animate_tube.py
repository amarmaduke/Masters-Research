import sys
import json
import math
import time
import numpy as np
from matplotlib import collections, transforms
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.path as mpath
import matplotlib.animation as animation

mpl.rcParams['axes.linewidth'] = .75
mpl.rcParams['text.usetex'] = True

def gen_lines(count, data) :
    n = data['n']
    m = data['m']
    S = data['sub_count'] + 1

    F = n*m + 1
    lines = np.zeros((F, 2, count))
    slines = np.zeros((S, 2, count))
    for c in range(count) :
        index = 'sindex%dtq%d' % (0, c)
        stuff = data[index]
        a = np.asarray(stuff)

        k, k2 = 0, 0
        for j in range(m) :
            x2, y2 = delta[j], 0
            
            lines[k][0][c], lines[k][1][c] = x2, y2
            k = k + 1
            for i in range(n) :
                x1, y1 = a[i + j*n], a[i + j*n + N]
                lines[k][0][c], lines[k][1][c] = x1, y1
                k = k + 1
                x2, y2 = x1, y1

            top_substrate_step = data['sub_h']
            top_substrate_count = data['sub_count']
            #bottom_substrate_step = data['osub_h']
            #bottom_substrate_count = data['osub_count']
            #bottom_substrate_x = data['osub']

            x2, y2 = a[2*N], a[2*N + 1]
            slines[k2][0][c], slines[k2][1][c] = x2, y2
            k2 = k2 + 1
            for j in range(top_substrate_count) :
                x1, y1 = x2 + top_substrate_step, y2
                slines[k2][0][c], slines[k2][1][c] = x1, y1
                x2, y2 = x1, y1
                k2 = k2 + 1

            #x2, y2 = bottom_substrate_x, 0
            #for j in range(bottom_substrate_count) :
            #    x1, y1 = x2 + bottom_substrate_step, y2
            #    plt.plot([x1, x2], [y1, y2], color='k', marker='o', markersize=marker_size)
            #    x2, y2 = x1, y1
    return lines, slines

json_raw = sys.stdin.read()
data = json.loads(json_raw)
n = data['n']
m = data['m']
print(m, n)
delta = data['delta']
sindex_count = data['sindex_count']
subh = data['sub_h']
S = data['sub_count']
N = n*m

tstep = data['save_step']
tmax = sindex_count*tstep
marker_size = 2

fig = plt.figure()

lines, slines = gen_lines(sindex_count, data)

l, s, = plt.plot([], [], 'ok-', [], [], 'ok-')
plt.xlim((-100, 0))
plt.ylim((0, 2))

def init() :
    l.set_data([], [])
    s.set_data([], [])
    return l, s,

def animate(i) :
    thisx = lines[0:96, 0, i]
    thisy = lines[0:96, 1, i]
    sthisx = slines[0:S, 0, i]
    sthisy = slines[0:S, 1, i]

    l.set_data(thisx, thisy)
    s.set_data(sthisx, sthisy)
    return l, s,

ani = animation.FuncAnimation(fig, animate, 1000, interval=25, blit=True, init_func=init)

plt.show()