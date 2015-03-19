import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = .75
mpl.rcParams['text.usetex'] = True

json_raw = sys.stdin.read()
data = json.loads(json_raw)
n = data['n']
m = data['m']
print(m, n)
delta = data['delta']
N = n*m

stuff = data['grid']
grid = np.asarray(stuff)
grid = grid[grid[:,0].argsort()]

K, K1 = np.shape(grid)
print(K)

reflect = 0

plt.figure(1, (3*8/4, 3*6/4))

mag_threshold = 10000000
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

    if k != 0 and magnitude < mag_threshold :
      mag_threshold = magnitude

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
        plt.plot(theta,mag,marker='*',markersize=5,color='k')


ax = plt.gca()
plt.scatter(X,Y,12,C,'o','Greys', vmin=0, vmax=96, linewidths=.3)

plt.xticks([0, 45, 90, 135, 180], fontsize=10)
plt.xlabel('$\\theta$', fontsize=12)
plt.yticks(fontsize=10)
plt.ylabel('$\sqrt{\lambda^2 + \mu^2}$', fontsize=12)
plt.axis('tight')
cbar = plt.colorbar(ticks=[0, 16, 32, 48, 64, 80, 96], fraction=.1)
cbar.solids.set_edgecolor('face')

cbar.ax.tick_params(labelsize=10)
ycticks = cbar.ax.yaxis.get_major_ticks()
ycticks[0].tick2On = False
ycticks[-1].tick2On = False
plt.xlim((0, 180))

xticks = ax.xaxis.get_major_ticks()
xticks[0].tick1On = False
xticks[0].tick2On = False
xticks[-1].tick1On = False
xticks[-1].tick2On = False

ax = plt.gca()
ax.invert_xaxis()

plt.savefig('temp_pull.eps', format='eps', dpi=300)
plt.show()
