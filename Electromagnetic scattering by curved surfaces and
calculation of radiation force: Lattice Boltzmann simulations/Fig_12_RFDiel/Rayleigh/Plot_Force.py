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


ratio = np.arange(0.02, 0.1+0.0001, 0.001)
ratioLBM = np.arange(0.02, 0.1+0.01, 0.01)



FxEr2LBM  = np.loadtxt('FxLBM_er_2.txt')
FxEr4LBM  = np.loadtxt('FxLBM_er_4.txt')
FxEr5LBM  = np.loadtxt('FxLBM_er_5.txt')



FxEr2Exact  = np.loadtxt('FxExact_er_2.txt')
FxEr4Exact  = np.loadtxt('FxExact_er_4.txt')
FxEr5Exact  = np.loadtxt('FxExact_er_5.txt')

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



ax.plot(ratio, FxEr2Exact, 'k-')
ax.plot(ratio, FxEr4Exact, 'r-')
ax.plot(ratio, FxEr5Exact, 'b-')


ax.plot(ratioLBM, FxEr2LBM, 'ko')
ax.plot(ratioLBM, FxEr4LBM, 'ro')
ax.plot(ratioLBM, FxEr5LBM, 'bo')


er2  = patches.Patch(color='k', label=r'$\varepsilon_r = 2$')
er4  = patches.Patch(color='r', label=r'$\varepsilon_r = 4$')
er5  = patches.Patch(color='b', label=r'$\varepsilon_r = 5$')


legend = plt.legend(handles=[er2, er4, er5])
ax.add_artist(legend)
    


                    
ax.set_xlabel(r'$a / \lambda$')
ax.set_ylabel(r'$\frac{\left< F_x \right> /L}{\lambda \varepsilon_0 E_0^2}$')




plt.savefig('FxRayleigh.svg')
plt.close()
##################################################################################################

