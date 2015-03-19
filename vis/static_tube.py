import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from  matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

## Script Inputs
# sindex is for grid files, (t is not important, set it to 0)
# t is for single simulations with many steps (sindex should be 0)
sindex, t = 0, 2

# [x_min, x_max, y_min, y_max]
view = [-80, 20, -40, 40]

# Size of the marker, 2 is good for plots with several points
marker_size = 2

# Sub views have to be either identically 0 or an array of 4 numbers
sub_view1 = [-2, 10, -1, 4]
marker_size1 = 4
zoom1 = 5
zoom1_loc = 1
connector1_loc1 = 3
connector1_loc2 = 4

sub_view2 = [-45, -20, -1, 4]
marker_size2 = 4
zoom2 = 4
zoom2_loc = 4
connector2_loc1 = 1
connector2_loc2 = 2

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

################################################################################

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

fig, ax = plt.subplots()
ax.set_aspect('equal',adjustable='datalim')

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
        ax.plot([x1, x2], [y1, y2], 
            color='k', marker='o', markersize=marker_size)
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
    ax.plot([x1, x2], [y1, y2], color='k', marker='o', markersize=marker_size)
    if sub_view1 != 0 :
        axins.plot([x1, x2], [y1, y2], color='k', marker='o', markersize=marker_size1)
    if sub_view2 != 0 :
        axins2.plot([x1, x2], [y1, y2], color='k', marker='o', markersize=marker_size2)
    x2, y2 = x1, y1

x2, y2 = bottom_substrate_x, 0
for j in range(bottom_substrate_count) :
    x1, y1 = x2 + bottom_substrate_step, y2
    ax.plot([x1, x2], [y1, y2], color='k', marker='o', markersize=marker_size)
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
    plt.setp(axins, xticks=[], yticks=[])
    mark_inset(ax, axins, loc1=connector1_loc1, loc2=connector1_loc2, fc="none", ec="0.5")

if sub_view2 != 0 :
    x1, x2, y1, y2 = sub_view2
    axins2.set_xlim(x1, x2)
    axins2.set_ylim(y1, y2)
    plt.setp(axins2, xticks=[], yticks=[])
    mark_inset(ax, axins2, loc1=connector2_loc1, loc2=connector2_loc2, fc="none", ec="0.5")

#plt.axis('equal')
#plt.axis(view)
plt.draw()
plt.show()
