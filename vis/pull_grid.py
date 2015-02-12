import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

## Script Inputs

mag_threshold = 200

################################################################################

def make_grayscale_colorlist_inc(N):
    color_list = []
    init_val = 0
    for i in range(N):
        color_list.append([init_val, init_val, init_val])
        init_val = min(init_val + 1/(N-1),1)
    return color_list

def make_grayscale_colorlist_dec(N):
    color_list = []
    init_val = 1
    for i in range(N):
        color_list.append([init_val, init_val, init_val])
        init_val = max(init_val - 1/(N-1),0)
    return color_list

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

gray_rolled = make_grayscale_colorlist_dec(N)
#gray_rolled = np.fliplr(gray_rolled)
my_cmap = colors.ListedColormap(gray_rolled, name="user_gray",N=N)

for k in range(K):
    i = grid[k,:]
    theta = i[0]
    mag = np.sqrt(i[1]**2 + i[2]**2)
    #colorList(i[7]+1)
    if i[3] == 1 and mag > mag_threshold:
        plt.plot(theta,mag,marker='o',markersize=3,color=gray_rolled[int(i[6])-1])
    elif i[3] == 2:
        plt.plot(theta,mag,marker='*',markersize=4,color='k')

ax = plt.gca()
ax.invert_xaxis()

# Much hack, wow. Fix with Line Collection?
sm = plt.cm.ScalarMappable(cmap=my_cmap, norm=plt.normalize(vmin=0, vmax=1))
sm._A = []
cb = plt.colorbar(sm, ticks=[0, .33, .66, 1])
cb.ax.set_yticklabels(['0', '32', '64', '96'])

plt.show()
