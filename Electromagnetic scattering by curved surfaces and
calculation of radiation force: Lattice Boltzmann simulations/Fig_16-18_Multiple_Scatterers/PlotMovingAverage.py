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

U_PEC_LD  = np.loadtxt('energy_ratio_1_a_40_PEC_LD.txt')
U_er5_LD  = np.loadtxt('energy_ratio_1_a_40_er5_LD.txt')



maEnergy_PEC = []
maEnergy_er5 = []

maEnergy_PEC_LD = []
maEnergy_er5_LD = []




step = 1

U = 1600




i = 0

while (i < len(U_PEC) - period + 1):
    window = U_PEC[i : i + period]
    window_average = sum(window) / period
    maEnergy_PEC.append(window_average / U)
    i += step

xmaEnergy_PEC = np.arange(0, len(maEnergy_PEC))/period

#################################


i = 0

while (i < len(U_er5) - period + 1):
    window = U_er5[i : i + period]
    window_average = sum(window) / period
    maEnergy_er5.append(window_average / U)
    i += step

xmaEnergy_er5 = np.arange(0, len(maEnergy_er5))/period

###################################




i = 0

while (i < len(U_PEC_LD) - period + 1):
    window = U_PEC_LD[i : i + period]
    window_average = sum(window) / period
    maEnergy_PEC_LD.append(window_average / U)
    i += step

xmaEnergy_PEC_LD = np.arange(0, len(maEnergy_PEC_LD))/period

#################################


i = 0

while (i < len(U_er5_LD) - period + 1):
    window = U_er5_LD[i : i + period]
    window_average = sum(window) / period
    maEnergy_er5_LD.append(window_average / U)
    i += step

xmaEnergy_er5_LD = np.arange(0, len(maEnergy_er5_LD))/period

###################################



xEnergy_PEC = np.arange(0, len(U_PEC))/period
xEnergy_er5 = np.arange(0, len(U_er5))/period

xEnergy_PEC_LD = np.arange(0, len(U_PEC_LD))/period
xEnergy_er5_LD = np.arange(0, len(U_er5_LD))/period


print(len(maEnergy_PEC))
print(np.absolute((maEnergy_PEC[1000*period] - maEnergy_PEC[400*period]) / maEnergy_PEC[400*period])*100)
print(np.absolute((maEnergy_er5[1000*period] - maEnergy_er5[400*period]) / maEnergy_er5[400*period])*100)


solpe_U_PEC = []
solpe_U_er5 = []

for i in range(200*period, 1000*period, period):

    solpe_U_PEC.append((maEnergy_PEC[i] - maEnergy_PEC[i-1]) / period)
    solpe_U_er5.append((maEnergy_er5[i] - maEnergy_er5[i-1]) / period)



##############################################################################################################################################################################


plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)





##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2.15), dpi=600, constrained_layout = True)


##ax1 = ax.inset_axes([0.4, 0.55, 0.55, 0.35])



ax.plot(xmaEnergy_PEC, maEnergy_PEC, 'k-')
ax.plot(xmaEnergy_er5, maEnergy_er5, 'r-')
##
##ax.plot(xmaEnergy_PEC_LD, maEnergy_PEC_LD, 'k--')
##ax.plot(xmaEnergy_er5_LD, maEnergy_er5_LD, 'r--')

##ax1.plot(np.arange(len(solpe_U_PEC)), solpe_U_PEC, 'k-')
##ax1.plot(np.arange(len(solpe_U_er5)), solpe_U_er5, 'r--')
##
##ax1.set_xticks(np.linspace(0, 800, 3))
##ax1.set_xticklabels([200, 600, 1000])
##
##ax1.set_xlabel(r'$\hat{t}$')
##ax1.set_ylabel(r'$\dot{\mathcal{U}}$')

##ax1.set_xscale('log')
##ax1.set_yscale('log')

##ax1.set_ylim(-0.001, 0.001)


ax.set_xlabel(r'$t / T$')
##ax.set_ylabel(r'$\mathcal{U}$')
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

plt.xscale('log')
##plt.grid()

plt.savefig('energy_multi_particles.svg')
plt.close(fig)

##################################################################################################
