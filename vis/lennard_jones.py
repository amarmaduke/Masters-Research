import sys
import json
import numpy as np
from matplotlib import collections, transforms
import matplotlib.pyplot as plt

eps = 100
sigma = 1

def lj_energy(x, eps, sigma):
    return eps*((sigma/x)**12 - 2*(sigma/x)**6)

def lj_force(x, eps, sigma):
    return (12*eps/sigma)*((sigma/x)**13 - (sigma/x)**7)

def lj_cut_energy(x, eps, sigma, cut):
    out = [] 
    for i in range(len(x)):
        if x[i] <= cut:
            temp = ( lj_energy(x[i], eps, sigma) 
                    - lj_energy(cut, eps, sigma)
                    + lj_force(cut, eps, sigma)*(x[i] - cut) )
            out.append(temp)
        else:
            out.append(0)
    return out

def lj_cut_force(x, eps, sigma, cut):
    out = [] 
    for i in range(len(x)):
        if x[i] <= cut:
            temp = ( lj_force(x[i], eps, sigma)
                    - lj_force(cut, eps, sigma) )    
            out.append(temp)
        else:
            out.append(0)
    return out



x = np.linspace(.8,10,10000)
reg_energy = lj_energy(x, eps, sigma)
reg_force = lj_force(x, eps, sigma)

rc = [1.5, 2, 2.5]
cut_energy = []
cut_force = []

for i in range(len(rc)):
   temp_energy = lj_cut_energy(x, eps, sigma, rc[i])
   temp_force = lj_cut_force(x, eps, sigma, rc[i])
   cut_energy.append(temp_energy)
   cut_force.append(temp_force)

plt.figure(1, (6, 5))
plt.axis([0, 3, -110, 110])
plt.plot(x,reg_energy, 'k')
plt.xticks([0, 1, 2, 3], fontsize=18)
plt.yticks([-100, 0, 100], fontsize=18)
plt.vlines(1, -110, 110, linestyles='dotted')
plt.annotate('$\sigma$', (1.05, 90), fontsize=28)
plt.hlines(-100, 0, 3, linestyles='dotted')
plt.annotate('$\\varepsilon$', (.1, -95), fontsize=28)
#plt.ylabel('$U(x)$', fontsize=24)

plt.savefig('lj1.eps', format='eps', dpi=1200)

plt.figure(2, (6, 5))
plt.axis([0, 3, -350, 350])
plt.plot(x,reg_force, 'k')
plt.xticks([0, 1, 2, 3], fontsize=18)
plt.yticks([-300, -200, -100, 0, 100, 200, 300], fontsize=18)
#plt.ylabel('-$U\'(x)$', fontsize=24)

plt.savefig('lj2.eps', format='eps', dpi=1200)

plt.figure(3, (6, 5))
plt.axis([0, 3, -110, 110])
plt.xticks([0, 1, 2, 3], fontsize=18)
plt.yticks([-100, 0, 100], fontsize=18)
plt.vlines(1, -110, 110, linestyles='dotted')
plt.annotate('$\sigma$', (1.05, 90), fontsize=28)
plt.hlines(-100, 0, 3, linestyles='dotted')
plt.annotate('$\\varepsilon$', (.1, -95), fontsize=28)

alpha_step = (.8 - .05)/len(rc)
for i in range(len(rc)):
    a = 1 - min(.2 + i*alpha_step,1)
    plt.plot(x, cut_energy[i], color=(a,a,a), label='$r_{cut} = %.1f$'%rc[i])
plt.plot(x, reg_energy, 'k')
plt.legend(fontsize=22)

plt.savefig('lj3.eps', format='eps', dpi=1200)

plt.figure(4, (6, 5))
plt.axis([0, 3, -350, 350])
plt.plot(x,reg_force, 'k')
plt.xticks([0, 1, 2, 3], fontsize=18)
plt.yticks([-300, -200, -100, 0, 100, 200, 300], fontsize=18)

alpha_step = (.8 - .05)/len(rc)
for i in range(len(rc)):
    a = 1 - min(.2 + i*alpha_step,1)
    plt.plot(x, cut_force[i], color=(a,a,a), label='$r_{cut} = %.1f$'%rc[i])
plt.plot(x, reg_force, 'k')
plt.legend(fontsize=22)

plt.savefig('lj4.eps', format='eps', dpi=1200)
plt.show()