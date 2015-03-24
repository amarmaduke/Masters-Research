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
    bS = data['osub_count'] + 1

    F = n*m + 1
    lines = np.zeros((F, 2, count))
    circles = np.zeros((F))
    slines = np.zeros((S, 2, count))
    blines = np.zeros((bS, 2))
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
            bottom_substrate_step = data['osub_h']
            bottom_substrate_count = data['osub_count']
            bottom_substrate_x = data['osub']

            x2, y2 = a[2*N], a[2*N + 1]
            slines[k2][0][c], slines[k2][1][c] = x2, y2
            k2 = k2 + 1
            for j in range(top_substrate_count) :
                x1, y1 = x2 + top_substrate_step, y2
                slines[k2][0][c], slines[k2][1][c] = x1, y1
                x2, y2 = x1, y1
                k2 = k2 + 1

            k3 = 0
            x2, y2 = bottom_substrate_x, 0
            blines[k3][0], blines[k3][1] = x2, y2
            k3 = k3 + 1
            for j in range(bottom_substrate_count) :
                x1, y1 = x2 + bottom_substrate_step, y2
                blines[k3][0], blines[k3][1] = x1, y1
                x2, y2 = x1, y1
                k3 = k3 + 1
    return lines, slines, blines

json_raw = sys.stdin.read()
data = json.loads(json_raw)
n = data['n']
m = data['m']
print(m, n)
delta = data['delta']
sindex_count = data['sindex_count']
subh = data['sub_h']
S = data['sub_count']
bS = data['osub_count']
N = n*m

tstep = data['save_step']
tmax = sindex_count*tstep
marker_size = 2

fig = plt.figure()

lines, slines, blines = gen_lines(sindex_count, data)

l, s, = plt.plot([], [], 'ok-', [], [], 'ok-')
plt.xlim((-100, 100))
plt.ylim((0, 100))

def init() :
    l.set_data([], [])
    s.set_data([], [])
    return l, s,

def animate(i) :
    print(i)
    thisx = lines[0:96, 0, i]
    thisy = lines[0:96, 1, i]
    sthisx = slines[0:S, 0, i]
    sthisy = slines[0:S, 1, i]

    l.set_data(thisx, thisy)
    s.set_data(sthisx, sthisy)
    return l, s,

#ani = animation.FuncAnimation(fig, animate, sindex_count, interval=1, blit=True, init_func=init)

'''
    The animation plots are made with the small scripts below and manually
    manuvered.
'''
g1 = (.6, .6, .6)
g2 = (.4, .4, .4)
g3 = (0, 0, 0)

#'''
# eb0.1/seed_t72_M38.5
fig.set_size_inches(3.25, 2.75,forward=True)
v = [-13, 1, -1, 15]
plt.axis(v)

plt.plot(lines[0:96, 0, 6000], lines[0:96, 1, 6000], 'o-', color=g1, markeredgecolor=g1)
plt.plot(slines[0:S, 0, 6000], slines[0:S, 1, 6000], 'o-', color=g1, markeredgecolor=g1)


plt.plot(lines[0:96, 0, 7000], lines[0:96, 1, 7000], 'o-', color=g2, markeredgecolor=g2)
plt.plot(slines[0:S, 0, 7000], slines[0:S, 1, 7000], 'o-', color=g2, markeredgecolor=g2)

plt.plot(lines[0:96, 0, 8000], lines[0:96, 1, 8000], 'o-', color=g3, markeredgecolor=g3)
plt.plot(slines[0:S, 0, 8000], slines[0:S, 1, 8000], 'o-', color=g3, markeredgecolor=g3)

plt.plot(blines[0:bS, 0], blines[0:bS, 1], 'o-', color=g3)

plt.xticks([-10, 0], fontsize=10)
plt.yticks([0, 10], fontsize=10)
plt.savefig('anim_temp.eps', format='eps', dpi=300)
#'''

'''
# eb0.1/seed_t160_M25.5
fig.set_size_inches(3.25, 2.75,forward=True)
v = [-60, 0, .75, 1.15]
plt.axis(v)

plt.plot(lines[0:96, 0, 6000], lines[0:96, 1, 6000], 'o-', color=g1, markeredgecolor=g1)
plt.plot(slines[0:S, 0, 6000], slines[0:S, 1, 6000], 'o-', color=g1, markeredgecolor=g1)


plt.plot(lines[0:96, 0, 7000], lines[0:96, 1, 7000], 'o-', color=g2, markeredgecolor=g2)
plt.plot(slines[0:S, 0, 7000], slines[0:S, 1, 7000], 'o-', color=g2, markeredgecolor=g2)

plt.plot(lines[0:96, 0, 8000], lines[0:96, 1, 8000], 'o-', color=g3, markeredgecolor=g3)
plt.plot(slines[0:S, 0, 8000], slines[0:S, 1, 8000], 'o-', color=g3, markeredgecolor=g3)

plt.plot(blines[0:bS, 0], blines[0:bS, 1], 'ok-', alpha=1)


plt.xticks([-50, -10], fontsize=10)
plt.yticks([.8, 1], fontsize=10)
plt.savefig('anim_temp.eps', format='eps', dpi=300)
#'''


plt.show()
