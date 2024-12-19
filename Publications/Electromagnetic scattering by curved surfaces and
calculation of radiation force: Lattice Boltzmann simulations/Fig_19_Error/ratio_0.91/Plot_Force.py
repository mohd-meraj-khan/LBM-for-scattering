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

ratioLBM = np.arange(0.02, 0.1 + 0.00001, 0.01)
ratioExact = np.arange(0.02, 0.1+0.000001, 0.00001)

############### loading LBM data ###################

fxExact = np.loadtxt(directory+'/FxExact.txt')
fxLBM = np.loadtxt(directory+'/FxLBM.txt')










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


plt.ion()


ax.plot(ratioExact, fxExact, 'k-', ratioLBM, fxLBM, 'ko')


ax.legend([r'Analytical', r'LBM'])
ax.set_xlabel(r'$a / \lambda$')
ax.set_ylabel(r'$\frac{\left< F_x \right>/L}{\lambda \varepsilon_0 E_0^2}$')

plt.savefig('Fx.svg')


plt.close()
##################################################################################################
