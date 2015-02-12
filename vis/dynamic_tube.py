import sys
import json
import numpy as np
from matplotlib import collections, transforms
import matplotlib.pyplot as plt

json_raw = sys.stdin.read()
data = json.loads(json_raw)
n = data['n']
m = data['m']
print(m, n)
delta = data['delta']
sindex_count = data['sindex_count']*.6
N = n*m

## Script Inputs
frames = 10

################################################################################

step = int(sindex_count / frames)

fig, axes = plt.subplots(1)

alpha_step = (.8 - .05)/frames

for i in range(frames):
    g = min(.05 + i*alpha_step,1)
    if i + 1 == frames:
        g = 1
    v = data['sindex0tq'+str(i*step)]
    x = delta + v[0:N]
    y = [0] + v[N:2*N]
    lines = list(zip(x, y))
    col = collections.LineCollection([lines], alpha=g, color='black')
    axes.add_collection(col, autolim=True)

axes.autoscale_view()
plt.show()
