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


ratio = np.arange(0.02, 2.6+0.000001, 0.00001)




FxEr2Ray  = np.loadtxt('FxExactRayleigh_er_2.txt')
FxEr2Mie  = np.loadtxt('FxExactMie_er_2.txt')
FxEr2GOp  = np.loadtxt('FxExactGO_er_2.txt')

FxEr2 = np.concatenate((FxEr2Ray, FxEr2Mie, FxEr2GOp))


FxEr4Ray  = np.loadtxt('FxExactRayleigh_er_4.txt')
FxEr4Mie  = np.loadtxt('FxExactMie_er_4.txt')
FxEr4GOp  = np.loadtxt('FxExactGO_er_4.txt')

FxEr4 = np.concatenate((FxEr4Ray, FxEr4Mie, FxEr4GOp))


FxEr5Ray  = np.loadtxt('FxExactRayleigh_er_5.txt')
FxEr5Mie  = np.loadtxt('FxExactMie_er_5.txt')
FxEr5GOp  = np.loadtxt('FxExactGO_er_5.txt')

FxEr5 = np.concatenate((FxEr5Ray, FxEr5Mie, FxEr5GOp))




########################################################################################################################################################
########################################################################################################################################################
			
plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)



##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2), dpi=600, constrained_layout = True)

n = -1

ax.plot(ratio[:n], FxEr2[:n], 'k-')
ax.plot(ratio[:n], FxEr4[:n], 'r-')
ax.plot(ratio[:n], FxEr5[:n], 'b-')

plt.axvline(x = 0.1, color='k', linestyle='--')
plt.axvline(x = 2, color='k', linestyle='--')

plt.text(0.03, 1.75, r'Rayleigh regime', verticalalignment='center', rotation='vertical')

plt.text(0.35, 2.5, r'Mie regime', verticalalignment='center', rotation='horizontal')

plt.text(2.6, 0.65, r'GO regime', verticalalignment='center', rotation='vertical')


ax.set_xlabel(r'$a / \lambda$')
ax.set_ylabel(r'$\frac{\left< F_x \right>/L}{\lambda \varepsilon_0 E_0^2}$')



plt.xscale('log')


##plt.grid()

plt.savefig('FxExactLog.svg')
plt.close()
##################################################################################################

