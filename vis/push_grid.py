import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms

json_raw = sys.stdin.read()
data = json.loads(json_raw)
n = data['n']
m = data['m']
print(m, n)
delta = data['delta']
N = n*m

x_range = data['range']
y_range = data['range2']

shape = [int((x_range[2] - x_range[0])/x_range[1]),
         int((y_range[2] - y_range[0])/y_range[1])]

grid = np.asarray(data['grid'])

#grid = grid[grid[:,0].argsort()] # Wow such sort

X = grid[:,0]
Y = np.sqrt(grid[:,1]**2 + grid[:,2]**2)
Z = grid[:,6]

X = np.reshape(X, (shape[0], shape[1]))
Y = np.reshape(Y, (shape[0], shape[1]))
Z = np.reshape(Z, (shape[0], shape[1]))

levels = []
for i in range(97):
    levels.append(i)

#plt.pcolormesh(X, Y, Z, cmap="Greys")

plt.contourf(X, Y, Z, levels=levels, cmap="Greys")

#plt.imshow(Z, cmap="Greys", interpolation='gaussian', aspect='auto', extent=[4, 176, 2, 22], vmin=0, vmax=96, origin='lower')

plt.xticks([45, 90, 135], ('$\\frac{\pi}{4}$', '$\\frac{\pi}{2}$', '$\\frac{3\pi}{2}$'), fontsize=16)
plt.xlabel('$\\theta$', fontsize=20)
plt.yticks([2, 6, 10, 14, 18,  22])
plt.ylabel('$\sqrt{\lambda^2 + \mu^2}$', fontsize=20)
plt.axis('tight')
plt.colorbar()
plt.show()