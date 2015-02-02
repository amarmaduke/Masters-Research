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

grid = np.asarray(data['grid'])

X = grid[:,0]
Y = np.sqrt(grid[:,1]**2 + grid[:,2]**2)
Z = grid[:,5]

print(np.shape(X))

X = np.reshape(X, (100, 100))
Y = np.reshape(Y, (100, 100))
Z = np.reshape(Z, (100, 100))

#print(Z)


levels = []

for i in range(97) :
    levels.append(i)

#plt.pcolormesh(X, Y, Z, cmap="Greys")

plt.contourf(X, Y, Z, levels=levels, cmap="Greys")

plt.show()