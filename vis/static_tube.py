import sys
import json
import yaml
import numpy as np
import matplotlib.pyplot as plt
from  matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = .75
mpl.rcParams['text.usetex'] = True

stream = open(sys.argv[1], 'r')
params = yaml.load(stream)

## Script Inputs
# sindex is for grid files, (t is not important, set it to 0)
# t is for single simulations with many steps (sindex should be 0)
sindex = params['sindex']
t = params['t']

# [x_min, x_max, y_min, y_max]
view = params['view']
# Size of the marker, 2 is good for plots with several points
marker_size = params['marker_size']

# Sub views have to be either identically 0 or an array of 4 numbers
sub_view1 = params['sub_view1']
marker_size1 = params['marker_size1']
zoom1 = params['zoom1']
zoom1_loc = params['zoom1_loc']
connector1_loc1 = params['connector1_loc1']
connector1_loc2 = params['connector1_loc2']

sub_view2 = params['sub_view2']
marker_size2 = params['marker_size2']
zoom2 = params['zoom2']
zoom2_loc = params['zoom2_loc']
connector2_loc1 = params['connector2_loc1']
connector2_loc2 = params['connector2_loc2']

## Zoom Loc Observations
# 1 -> Top Right
# 2 -> Top Left
# 3 -> Bottom Left
# 4 -> Bottom Right
# 5 -> Center Right
# 6 -> Center Left
# 7 -> Right Center
# 8 -> Bottom Center
# 9 -> Top Center
# 10 -> Center
# >10 Error Code, <1 Error Code

json_raw = sys.stdin.read()
data = json.loads(json_raw)
n = data['n']
m = data['m']
print(m, n)
delta = data['delta']
N = n*m

index = 'sindex%dtq%d' % (sindex,t,)
stuff = data[index]
a = np.asarray(stuff)

plt.figure(1, (3.25, 2.75))
ax = plt.gca()

ax.set_aspect('equal',adjustable='datalim')

plt.xticks([0], fontsize=10)
plt.yticks([0], fontsize=10)

if sub_view1 != 0 :
    axins = zoomed_inset_axes(ax, zoom1, loc=zoom1_loc)
    axins.set_aspect('equal',adjustable='datalim')

if sub_view2 != 0 :
    axins2 = zoomed_inset_axes(ax, zoom2, loc=zoom2_loc)
    axins2.set_aspect('equal',adjustable='datalim')

for j in range(m) :
    x2, y2 = delta[j], 0
    for i in range(n) :
        x1, y1 = a[i + j*n], a[i + j*n + N]
        if marker_size != 0 :
            ax.plot([x1, x2], [y1, y2],
                color='k', marker='o', markersize=marker_size)
        else :
            ax.plot([x1, x2], [y1, y2], color='k')
        if sub_view1 != 0 :
            axins.plot([x1, x2], [y1, y2],
                color='k', marker='o', markersize=marker_size1)
        if sub_view2 != 0 :
            axins2.plot([x1, x2], [y1, y2],
                color='k', marker='o', markersize=marker_size2)
        x2, y2 = x1, y1

top_substrate_step = data['sub_h']
top_substrate_count = data['sub_count']
bottom_substrate_step = data['osub_h']
bottom_substrate_count = data['osub_count']
bottom_substrate_x = data['osub']

x2, y2 = a[2*N], a[2*N + 1]
for j in range(top_substrate_count) :
    x1, y1 = x2 + top_substrate_step, y2
    if marker_size != 0 :
        ax.plot([x1, x2], [y1, y2], color='k', marker='o', markersize=marker_size)
    else :
        ax.plot([x1, x2], [y1, y2], color='k')
    if sub_view1 != 0 :
        axins.plot([x1, x2], [y1, y2], color='k', marker='o', markersize=marker_size1)
    if sub_view2 != 0 :
        axins2.plot([x1, x2], [y1, y2], color='k', marker='o', markersize=marker_size2)
    x2, y2 = x1, y1

x2, y2 = bottom_substrate_x, 0
for j in range(bottom_substrate_count) :
    x1, y1 = x2 + bottom_substrate_step, y2
    if marker_size != 0 :
        ax.plot([x1, x2], [y1, y2], color='k', marker='o', markersize=marker_size)
    else :
        ax.plot([x1, x2], [y1, y2], color='k')
    if sub_view1 != 0 :
        axins.plot([x1, x2], [y1, y2], color='k', marker='o', markersize=marker_size1)
    if sub_view2 != 0 :
        axins2.plot([x1, x2], [y1, y2], color='k', marker='o', markersize=marker_size2)
    x2, y2 = x1, y1

x1, x2, y1, y2 = view
ax.set_xlim(x1, x2)
ax.set_ylim(y1, y2)

if sub_view1 != 0 :
    x1, x2, y1, y2 = sub_view1
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    plt.setp(axins, xticks=[0], yticks=[0])
    axins.set_xticklabels(' ', fontsize=1)
    axins.set_yticklabels(' ', fontsize=1)
    mark_inset(ax, axins, loc1=connector1_loc1, loc2=connector1_loc2, fc="none", ec="0.75")

if sub_view2 != 0 :
    x1, x2, y1, y2 = sub_view2
    axins2.set_xlim(x1, x2)
    axins2.set_ylim(y1, y2)
    plt.setp(axins2, xticks=[0], yticks=[0])
    axins2.set_xticklabels(' ', fontsize=1)
    axins2.set_yticklabels(' ', fontsize=1)
    mark_inset(ax, axins2, loc1=connector2_loc1, loc2=connector2_loc2, fc="none", ec="0.75")

plt.draw()

plt.savefig('static_temp.eps', format='eps', dpi=300)
plt.show()
