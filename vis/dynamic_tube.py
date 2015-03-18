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
N = n*m

tstep = data['save_step']
tmax = data['tmax']

## Script Inputs
frames = 500
frame_begin = 1
frame_stop = 10000
view = [-96, 20, -2, 2]

particle_begin = 0
particle_end = 16
particle_step = 1

#particle_begin = 10
#particle_end = 19
#particle_step = 1

#particle_begin = 96-18
#particle_end = 96-9
#particle_step = 1

#particle_begin = 96-9
#particle_end = 96
#particle_step = 1

################################################################################

stop = min(sindex_count, frame_stop - frame_begin)
step = int(stop / frames)

print(sindex_count)
linestyles = ['-', '--', '-.', ':']

alpha_step = (1 - 0.05)/frames

counter = 0
#for i in range(frames):
#    g = .05 + i*alpha_step
#    if i + 1 == frames:
#        g = 1
#    v = data['sindex0tq'+str(frame_begin + i*step)]
    #x = delta + v[0:N-1]
#    x = v[particle_begin:particle_end:particle_step] #+ v[0:N-1:int((N-2)/8)]
    #y = [0] + v[N:2*N-1]
#    y = v[N+particle_begin:N+particle_end:particle_step] #+ v[N:2*N-1:int((N-2)/8)]
#    c = [0*i + 5 for i in range(1, N+1)]
#    lines = list(zip(x, y))
#    col = collections.LineCollection([lines], alpha=g, color='black')
    #plt.scatter(x, y, s=3, alpha=g, color='black')
#    axes.add_collection(col, autolim=True)
#    counter = counter + 1

sub_count = int(math.ceil(math.sqrt(particle_end - particle_begin)/particle_step))

f, axarr = plt.subplots(sub_count, sub_count, squeeze=True, figsize=(6,5))
ii = 0
jj = 0
for i in range(particle_begin, particle_end, particle_step) :
    g = 1 - counter/(particle_end - particle_begin)
    t = [j*tstep for j in range(0, sindex_count)]
    x = []
    y = []
    for j in range(0, sindex_count) :
        v = data['sindex0tq'+str(j)]
        x.extend([v[i]])
        y.extend([v[i+N]])
    plt.figure(2)
    plt.plot(t, y, alpha=g, color='black', linestyle=linestyles[counter % 4])
    plt.plot(t[0], y[0], alpha=g, color='black', marker='D', markersize=8)
    plt.axis('tight')
    plt.figure(3)
    plt.plot(t, x, alpha=g, color='black', linestyle=linestyles[counter % 4])
    plt.plot(t[0], x[0], alpha=g, color='black', marker='D', markersize=8)
    plt.axis('tight')
    #plt.figure(4)
    #plt.plot(x, y, alpha=g, color='black', linestyle=linestyles[counter % 4])
    #plt.plot(x[0], y[0], alpha=g, color='black', marker='D', markersize=8)
    #plt.axis('tight')

    axarr[ii, jj].plot(x, y, alpha=1, linestyle=linestyles[0], color='black')
    #axarr[ii, jj].add_collection(col)

    axarr[ii, jj].plot(x[0], y[0], alpha=1, color='black', marker='D', markersize=4)
    axarr[ii, jj].axis('tight')
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
    xt = [xt[-3], xt[2]]
    yt = axarr[ii, jj].get_yticks()
    yt = [yt[-3], yt[2]]
    plt.setp(axarr[ii, jj], xticks=xt, yticks=yt)
    axarr[ii, jj].get_xaxis().get_major_formatter().set_useOffset(False)
    axarr[ii, jj].get_yaxis().get_major_formatter().set_useOffset(False)
    axarr[ii, jj].tick_params(axis='both', which='major', labelsize=10)

    jj = jj + 1
    if jj >= sub_count :
        ii = ii + 1
        jj = 0
    counter = counter + 1


plt.tight_layout()

plt.show()
