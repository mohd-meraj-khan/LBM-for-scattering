import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.patches as patches
import sys
import os
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset


directory = 'data'

ratioMie = np.arange(0.5, 4 + 0.01, 0.01)
ratioRayleigh = np.arange(0.02, 0.1 + 0.0001, 0.001)

############### loading LBM data ###################

fxRayleighExact = np.loadtxt(directory+'/FxRayleighExact.txt')
fxMieExact = np.loadtxt(directory+'/FxMieExact.txt')



fxLBMMie1 = np.loadtxt(directory+'/FxAvgMie.txt')
fxLBMMie2 = np.loadtxt(directory+'/FxAvgGO.txt')

fxLBMMie = np.concatenate((fxLBMMie1, fxLBMMie2))

fxLBMRayleigh = np.loadtxt(directory+'/FxAvgRayleigh.txt')



ratioLBMMie = np.arange(0.5, 4 + 0.01, 0.5)
ratioLBMRayleigh = np.arange(0.02, 0.1 + 0.001, 0.01)


############################################################################################################################################################################

c = 2

plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)



##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2), dpi=600, constrained_layout = True)

ax1 = ax.inset_axes([0.65, 0.15, 0.32, 0.35])

plt.ion()


ax.plot(ratioMie, fxMieExact, 'k-', ratioLBMMie, fxLBMMie, 'ko')
ax1.plot(ratioRayleigh, fxRayleighExact, 'k-', ratioLBMRayleigh, fxLBMRayleigh, 'ko')

ax.legend([r'Analytical', r'LBM'])
ax.set_xlabel(r'$a / \lambda$')
ax.set_ylabel(r'$\frac{\left< F_x \right>/L}{\lambda \varepsilon_0 E_0^2}$')

plt.savefig('FxPEC.svg')


plt.close()
##################################################################################################
