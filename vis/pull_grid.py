import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

## Script Inputs

mag_threshold = 0

################################################################################

json_raw = sys.stdin.read()
data = json.loads(json_raw)
n = data['n']
m = data['m']
print(m, n)
delta = data['delta']
N = n*m

stuff = data['grid']
grid = np.asarray(stuff)

grid = grid[grid[:,0].argsort()] # Wow such sort

K, K1 = np.shape(grid)
print(K)

reflect = 0

theta_p = -1
theta_c = 0
counter = 0
M = 0
m = 1000000
M_adhered = []
m_pulloff = []
theta_v = []
for k in range(K) :
    i = grid[k,:]

    magnitude = np.sqrt(i[1]**2 + i[2]**2);
    theta_c = i[0];

    if k == 0 :
        theta_p = i[0]
    
    if theta_c != theta_p or k == K-1 :
        M_adhered.insert(counter, M)
        m_pulloff.insert(counter, m)
        #if reflect and theta_p > 90 :
        #    theta_v.insert(counter, 180 - theta_p)
        #else :
        theta_v.insert(counter, theta_p)
        M = 0
        m = 1000000
        counter = counter + 1;
        theta_p = theta_c
    
    if i[3] == 1 :
        M = max(magnitude,M)
    elif i[3] == 0 :
        m = min(magnitude,m)

plt.plot(theta_v,M_adhered,'--',color='k')
plt.plot(theta_v,m_pulloff,color='k')
plt.plot(theta_v,m_pulloff,color='k')

X = []
Y = []
C = []

for k in range(K):
    i = grid[k,:]
    theta = i[0]
    mag = np.sqrt(i[1]**2 + i[2]**2)
    #colorList(i[7]+1)
    if i[3] == 1 and mag > mag_threshold:
        X.append(theta);
        Y.append(mag);
        C.append(i[6])
    elif i[3] == 2:
        plt.plot(theta,mag,marker='*',markersize=7,color='k')

plt.scatter(X,Y,12,C,'o','Greys', vmin=0, vmax=96, linewidths=.25)

plt.colorbar()
plt.axis('tight')
plt.xlim((0, 180))

ax = plt.gca()
ax.invert_xaxis()

plt.savefig('temp_pull.eps', format='eps', dpi=1200)
plt.show()
