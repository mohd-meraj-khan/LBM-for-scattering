import numpy as np ; from numpy.linalg import *
import sys
import os
import csv
import math
import scipy.special as sc
from math import e
import cmath
import matplotlib.pyplot as plt
import matplotlib.patches as patches


ratio = 2

a = np.arange(10, 100+1, 10)


step = 1

FxLBM  = np.loadtxt('data/FxLBM.txt')

U = []
U_MA = []
xU_MA = []

for i in range(len(a)):
    u = np.loadtxt('data/energy_ratio_2_a_{}.txt'.format(a[i]))
    U.append(u / a[i]**2)

    period = int(np.round(3 * a[i] / ratio))

    uma = []

    j = 0
    while (j < len(u) - period + 1):
        window = u[j : j + period]
        window_average = sum(window) / period
        uma.append(window_average / a[i]**2)
        j += step
        
    U_MA.append(uma)
    xU_MA.append(np.arange(0, len(uma)))


########################################################################################################################################################
########################################################################################################################################################
			
plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 1)
plt.rc('text', usetex = True)



##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2), dpi=600, constrained_layout = True)

ax.plot(a/ratio, FxLBM[10:], 'k-o')
                 
ax.set_xlabel(r'$\lambda / \Delta x$')
ax.set_ylabel(r'$\frac{\left< F_x \right> /L}{\lambda \varepsilon_0 E_0^2}$')


plt.savefig('FxHexagon.svg')
plt.close()
##################################################################################################




##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2.25), dpi=600, constrained_layout = True)

for i in range(len(a)):
##    ax.plot(np.arange(len(U[i])), U[i])
    ax.plot(xU_MA[i], U_MA[i])
                 
ax.set_xlabel(r'$t / \Delta t$')
ax.set_ylabel(r'$\frac{\left< U \right> / L}{\varepsilon_0 E_0^2 a^2}$')

plt.legend([r'$\lambda / \Delta x = 5$', r'$\lambda / \Delta x = 10$', r'$\lambda / \Delta x = 15$', r'$\lambda / \Delta x = 20$',
            r'$\lambda / \Delta x = 25$', r'$\lambda / \Delta x = 30$', r'$\lambda / \Delta x = 35$', r'$\lambda / \Delta x = 40$',
            r'$\lambda / \Delta x = 45$', r'$\lambda / \Delta x = 50$'], ncol=2)

##plt.xscale('log')
##plt.yscale('log')

plt.savefig('UHexagon.svg')
plt.close()
##################################################################################################

