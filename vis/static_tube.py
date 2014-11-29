import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from  matplotlib.animation import FuncAnimation

json_raw = sys.stdin.read()
data = json.loads(json_raw)
n = data['n']
m = data['m']
print m, n
delta = data['delta']
N = n*m

index = 'sindex0tq%d' % (1000,)
stuff = data[index]
a = np.asarray(stuff);

for j in range(m) :
    x2, y2 = delta[j], 0
    for i in range(n) :
        x1, y1 = a[i + j*n], a[i + j*n + N]
        plt.plot([x1, x2], [y1, y2], color='k', marker='o', markersize=2)
        x2, y2 = x1, y1

top_substrate_step = data['sub_h']
top_substrate_count = data['sub_count']
bottom_substrate_step = data['osub_h']
bottom_substrate_count = data['osub_count']
bottom_substrate_x = data['osub']

x2, y2 = a[2*N], a[2*N + 1]
for j in range(top_substrate_count) :
    x1, y1 = x2 + top_substrate_step, y2
    plt.plot([x1, x2], [y1, y2], color='k', marker='o', markersize=2)
    x2, y2 = x1, y1

x2, y2 = bottom_substrate_x, 0
for j in range(bottom_substrate_count) :
    x1, y1 = x2 + bottom_substrate_step, y2
    plt.plot([x1, x2], [y1, y2], color='k', marker='o', markersize=2)
    x2, y2 = x1, y1



plt.axis('equal')
plt.show()