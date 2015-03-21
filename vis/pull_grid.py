import sys
import json
import math
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

mpl.rcParams['figure.figsize'] = 6, 3

fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.plot(theta_v, M_adhered, '--', color='k')
ax1.plot(theta_v, m_pulloff, color='k')

X = []
Y = []
C = []

for k in range(K):
    i = grid[k,:]
    theta = i[0]
    mag = np.sqrt(i[1]**2 + i[2]**2)
    if i[3] == 1 and mag > mag_threshold:
        X.append(theta);
        Y.append(mag);
        C.append(i[6])
    elif i[3] == 2:
        ax1.plot(theta,mag,marker='*',markersize=5,color='k')

ax1.scatter(X,Y,12,C,'o','Greys', vmin=0, vmax=96, linewidths=.3)

l_d = [-m_pulloff[k]*math.sin(math.pi*theta_v[k]/180) for k in range(len(theta_v))]
m_d = [m_pulloff[k]*math.cos(math.pi*theta_v[k]/180) for k in range(len(theta_v))]
l_a = [-M_adhered[k]*math.sin(math.pi*theta_v[k]/180) for k in range(len(theta_v))]
m_a = [M_adhered[k]*math.cos(math.pi*theta_v[k]/180) for k in range(len(theta_v))]
l_ = [-Y[k]*math.sin(math.pi*X[k]/180) for k in range(len(X))]
m_ = [Y[k]*math.cos(math.pi*X[k]/180) for k in range(len(X))]
 
ax2.plot(m_a, l_a, '--', color='k')
ax2.plot(m_d, l_d, color='k')
s1 = ax2.scatter(m_, l_, 12, C, 'o', 'Greys', vmin=0, vmax=96, linewidths=.3)
ax2.axis('tight')

ax1.invert_xaxis()
ax1.tick_params(axis='both', which='major', labelsize=10)
ax1.set_xticks([0, 45, 90, 135, 180])
ax1.set_xlabel('$\\theta$', fontsize=12)
ax1.set_ylabel('$\sqrt{\lambda^2 + \mu^2}$', fontsize=12)
ax1.axis('tight')
ax1.set_xlim((0, 180))

ax2.tick_params(axis='both', which='major', labelsize=10)
ax2.set_xlabel('$\mu$', fontsize=12)
ax2.set_ylabel('$\lambda$', fontsize=12)

plt.tight_layout()

xticks = ax1.xaxis.get_major_ticks()
xticks[0].tick1On = False
xticks[0].tick2On = False
xticks[-1].tick1On = False
xticks[-1].tick2On = False

cbar = fig.colorbar(s1, ax=[ax1, ax2], ticks=[0, 16, 32, 48, 64, 80, 96], fraction=.1)
cbar.solids.set_edgecolor('face')

cbar.ax.tick_params(labelsize=10)
ycticks = cbar.ax.yaxis.get_major_ticks()
ycticks[0].tick2On = False
ycticks[-1].tick2On = False

plt.savefig('temp_pull.eps', format='eps', dpi=300)
plt.show()
