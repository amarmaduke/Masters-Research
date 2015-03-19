import sys
import json
import math
import numpy as np
from matplotlib import collections, transforms
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.path as mpath

mpl.rcParams['axes.linewidth'] = .75
mpl.rcParams['text.usetex'] = True

json_raw = sys.stdin.read()
data = json.loads(json_raw)
n = data['n']
m = data['m']
print(m, n)
delta = data['delta']
sindex_count = data['sindex_count']
subh = data['sub_h']
N = n*m

tstep = data['save_step']
tmax = sindex_count*tstep
print(tmax)

## Script Inputs
frames = 500
frame_begin = 1
frame_stop = 10000
view = [-96, 20, -2, 2]

shape = [4, 3]

pbegin = 3
pstep = 1

################################################################################

pend = pbegin + shape[0]*shape[1]
stop = min(sindex_count, frame_stop - frame_begin)
step = int(stop / frames)

print(sindex_count)
linestyles = ['-', '--', '-.', ':']

alpha_step = (1 - 0.05)/frames

detach_t = 0
for i in range(0, sindex_count) :
    v = data['sindex0tq'+str(i)]
    attached = False
    for j in range(0, data['sub_count']) :
        for k in range(0, N) :
            sub_x = v[2*N]
            sub_y = v[2*N + 1]
            sx = sub_x + subh*j
            sy = sub_y
            x = v[k]
            y = v[k + N]
            distance = math.sqrt((x - sx)**2 + (y - sy)**2)
            tol = 1.123 + 1e-6
            if distance <= tol :
                attached = True
    if not attached :
        detach_t = i
        break


counter = 0

f, axarr = plt.subplots(shape[0], shape[1])
ii = 0
jj = 0
for i in range(pbegin, pend) :
    g = 1 - counter/(pend - pbegin)
    t = [j*tstep for j in range(0, sindex_count)]
    x = []
    y = []
    for j in range(0, sindex_count) :
        v = data['sindex0tq'+str(j)]
        x.extend([v[i*pstep]])
        y.extend([v[i*pstep + N]])

    axarr[ii, jj].plot(x, y, alpha=1, color='black')

    axarr[ii, jj].plot(x[0], y[0], color='black', marker='D', markersize=4)
    if (i + 1) % 10 == 1 :
        axarr[ii, jj].set_title('$%d$st particle'% (i+1), fontsize=12)
    elif (i + 1) % 10 == 2 :
        axarr[ii, jj].set_title('$%d$nd particle'% (i+1), fontsize=12)
    elif (i + 1) % 10 == 3 :
        axarr[ii, jj].set_title('$%d$rd particle'% (i+1), fontsize=12)
    else :
        axarr[ii, jj].set_title('$%d$th particle'% (i+1), fontsize=12)
    v = axarr[ii, jj].axis()
    xbuff = 0.05*abs(v[1] - v[0])
    ybuff = 0.05*abs(v[3] - v[2])
    axarr[ii, jj].axis([v[0]-xbuff, v[1]+xbuff, v[2]-ybuff, v[3]+ybuff])
    xt = axarr[ii, jj].get_xticks()
    xt = [xt[-2], xt[2]]
    yt = axarr[ii, jj].get_yticks()
    yt = [yt[-2], yt[2]]
    plt.setp(axarr[ii, jj], xticks=xt, yticks=yt)
    axarr[ii, jj].get_xaxis().get_major_formatter().set_useOffset(False)
    axarr[ii, jj].get_yaxis().get_major_formatter().set_useOffset(False)
    axarr[ii, jj].tick_params(axis='both', which='major', labelsize=10)

    jj = jj + 1
    if jj >= shape[1] :
        ii = ii + 1
        jj = 0
    counter = counter + 1

plt.tight_layout()

counter = 0
for i in range(pbegin, pend) :
    g = 1 - counter/(pend - pbegin)
    t = [j*tstep for j in range(0, sindex_count)]
    x = []
    y = []
    sx = []
    for j in range(0, sindex_count) :
        v = data['sindex0tq'+str(j)]
        x.extend([v[i*pstep]])
        y.extend([v[i*pstep + N]])
        sx.extend([v[2*N]])
    plt.figure(2)
    plt.title("$\Delta y by time$")
    plt.plot(t, y, alpha=g, color='black', linestyle=linestyles[counter % 4])
    plt.plot(t[0], y[0], alpha=g, color='black', marker='D', markersize=8)
    plt.plot(detach_t*tstep, y[detach_t], color='black', marker='o', markersize=8)
    plt.axis('tight')
    plt.figure(3)
    plt.title("$\Delta x (normalized) by time$")
    xx = [k + abs(x[0]) for k in x]
    plt.plot(t, xx, alpha=g, color='black', linestyle=linestyles[counter % 4])
    plt.plot(t[0], xx[0], alpha=g, color='black', marker='D', markersize=8)
    plt.plot(detach_t*tstep, xx[detach_t], color='black', marker='o', markersize=8)
    plt.axis('tight')
    plt.figure(4)
    plt.title("$\Delta sub_x by time$")
    plt.plot(t, sx, alpha=g, color='black', linestyle=linestyles[counter % 4])
    plt.axis('tight')
    counter = counter + 1

a, b, counter = 0, 0, 0
for i in range(0, detach_t, 5) :
#for i in range(int(5.0/tstep), detach_t) :
#for i in range(detach_t-100, detach_t, 1) :
    g = 1 - 4*counter/(detach_t)
    v = data['sindex0tq'+str(i)]
    attached = False
    d = []
    for k in range(1, N) :
        sx = v[k - 1]
        sy = v[k + N - 1]
        x = v[k]
        y = v[k + N]
        distance = math.sqrt((x - sx)**2 + (y - sy)**2)
        d.append(distance)

    a = b
    b = a + .0025
    #d = [k + a for k in d]
    t = [k for k in range(0, len(d))]
    plt.figure(5)
    plt.title("$\Delta pairwise distance by time$")
    plt.plot(d, t, color='black', alpha=g)
    counter = counter + 1

plt.show()
