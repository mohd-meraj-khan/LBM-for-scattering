import numpy as np ; from numpy.linalg import *
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.patches as patches
import sys
import os
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
import csv
import math
import scipy.special as sc
from math import e
import cmath
from scipy import interpolate
from matplotlib.gridspec import GridSpec

er = 5
a = 40





ratio = 1
period = int(np.round(3 * a / ratio))


U_PEC  = np.loadtxt('energy_ratio_1_a_40_PEC.txt')

U_er5  = np.loadtxt('energy_ratio_1_a_40_er5.txt')



maEnergy_PEC = []
maEnergy_er5 = []





step = 1

U = 1600




i = 0

while (i < len(U_PEC) - period + 1):
    window = U_PEC[i : i + period]
    window_average = sum(window) / period
    maEnergy_PEC.append(window_average / U)
    i += step

xmaEnergy_PEC = np.arange(0, len(maEnergy_PEC))/period

print(len(U_PEC))

#################################


i = 0

while (i < len(U_er5) - period + 1):
    window = U_er5[i : i + period]
    window_average = sum(window) / period
    maEnergy_er5.append(window_average / U)
    i += step

xmaEnergy_er5 = np.arange(0, len(maEnergy_er5))/period

###################################



xEnergy_PEC = np.arange(0, len(U_PEC))/period
xEnergy_er5 = np.arange(0, len(U_er5))/period

    


##############################################################################################################################################################################


plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)


n = 0*127

##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2.15), dpi=600, constrained_layout = True)


ax.plot(xmaEnergy_PEC[n:], maEnergy_PEC[n:], 'k-')
ax.plot(xmaEnergy_er5[n:], maEnergy_er5[n:], 'r-')


ax.set_xlabel(r'$t/T$')
ax.set_ylabel(r'$\frac{\left< U \right> / L}{\varepsilon_0 E_0^2 a^2}$')

plt.legend([r'PEC', r'$\varepsilon_r = 5$'])

##plt.xscale('log')

##plt.grid()

plt.savefig('MA_energy_multi_particles.svg')
plt.close(fig)

##################################################################################################



##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2.15), dpi=600, constrained_layout = True)


ax.plot(xEnergy_PEC, U_PEC/U, 'k-')
ax.plot(xEnergy_er5, U_er5/U, 'r-')



ax.set_xlabel(r'$t/T$')
ax.set_ylabel(r'$\frac{\left< U \right> / L}{\varepsilon_0 E_0^2 a^2}$')

plt.legend([r'PEC', r'$\varepsilon_r = 5$'])

##plt.xscale('log')
##plt.grid()

plt.savefig('energy_multi_particles.svg')
plt.close(fig)

##################################################################################################
