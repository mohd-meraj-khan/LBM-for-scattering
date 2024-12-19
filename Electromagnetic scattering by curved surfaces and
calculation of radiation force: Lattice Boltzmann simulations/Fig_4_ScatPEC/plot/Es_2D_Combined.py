import numpy as np ; from numpy.linalg import *
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import matplotlib.markers
import sys
import os
from matplotlib import patches

import matplotlib as mpl
from matplotlib import rc





#######################################################

Es_2D_005 = np.loadtxt('Ez_s_0.05.txt')
Es_2D_05 = np.loadtxt('Ez_s_0.5.txt')
Es_2D_1 = np.loadtxt('Ez_s_1.txt')
Es_2D_4 = np.loadtxt('Ez_s_4.txt')




c = 2

plt.clf()

plt.rc('font', family = 'serif', size = 8)
plt.rc('xtick', labelsize = 8)
plt.rc('ytick', labelsize = 8)
plt.rc('lines', markersize = 3, lw = 0.5)
plt.rc('text', usetex = True)


##################################################################################################




n005, n05, n1, n4 = 625, 625, 625, 1250










fig, axes = plt.subplots(nrows = 2, ncols = 2, figsize = (6.69/c, 7.4/c), constrained_layout=True, dpi=600)

im1 = axes[0,0].imshow(Es_2D_005, vmin = -1, vmax = 1, cmap='seismic', origin='lower')
im2 = axes[0,1].imshow(Es_2D_05, vmin = -1, vmax = 1, cmap='seismic', origin='lower')
im3 = axes[1,0].imshow(Es_2D_1, vmin = -1, vmax = 1, cmap='seismic', origin='lower')
im4 = axes[1,1].imshow(Es_2D_4, vmin = -1, vmax = 1, cmap='seismic', origin='lower')



axes[0,0].set_xticks(np.linspace(0, n005, 6))
axes[0,0].set_xticklabels([])
axes[0,0].set_yticks(np.linspace(0, n005, 6))
axes[0,0].set_yticklabels([0, 5, 10, 15, 20, 25])
rect1 = plt.Circle((n005/2, n005/2), 25, color='w', alpha=1, linewidth=1)
rect2 = plt.Circle((n005/2, n005/2), 25, color='k', alpha=0.5)
axes[0,0].add_patch(rect1)
axes[0,0].add_patch(rect2)
axes[0,0].set_ylabel(r'$y/a$', fontsize=10)
##axes[0,0].tick_params(which = 'both', direction='in')
axes[0,0].text(0.1*n005, 0.85*n005, r'$a/\lambda = 0.05$', backgroundcolor='white')

axes[0,1].set_xticks(np.linspace(0, n05, 6))
axes[0,1].set_xticklabels([])
axes[0,1].set_yticks(np.linspace(0, n05, 6))
axes[0,1].set_yticklabels([])
rect1 = plt.Circle((n05/2, n05/2), 25, color='w', alpha=1)
rect2 = plt.Circle((n05/2, n05/2), 25, color='k', alpha=0.5)
axes[0,1].add_patch(rect1)
axes[0,1].add_patch(rect2)
##axes[0,1].tick_params(which = 'both', direction='in')
axes[0,1].text(0.1*n05, 0.85*n05, r'$a/\lambda = 0.5$', backgroundcolor='white')

axes[1,0].set_xticks(np.linspace(0, n1, 6))
axes[1,0].set_xticklabels([0, 5, 10, 15, 20, 25])
axes[1,0].set_yticks(np.linspace(0, n1, 6))
axes[1,0].set_yticklabels([0, 5, 10, 15, 20, 25])
rect1 = plt.Circle((n1/2, n1/2), 25, color='w', alpha=1)
rect2 = plt.Circle((n1/2, n1/2), 25, color='k', alpha=0.5)
axes[1,0].add_patch(rect1)
axes[1,0].add_patch(rect2)
axes[1,0].set_xlabel(r'$x/a$', fontsize=10)
axes[1,0].set_ylabel(r'$y/a$', fontsize=10)
##axes[1,0].tick_params(which = 'both', direction='in')
axes[1,0].text(0.1*n1, 0.85*n1, r'$a/\lambda = 1$', backgroundcolor='white')

axes[1,1].set_xticks(np.linspace(0,n4,6))
axes[1,1].set_xticklabels([0, 5, 10, 15, 20, 25])
axes[1,1].set_yticks(np.linspace(0,n4,6))
axes[1,1].set_yticklabels([])
rect1 = plt.Circle((n4/2, n4/2), 50, color='w', alpha=1)
rect2 = plt.Circle((n4/2, n4/2), 50, color='k', alpha=0.5)
axes[1,1].add_patch(rect1)
axes[1,1].add_patch(rect2)
axes[1,1].set_xlabel(r'$x/a$', fontsize=10)
##axes[1,1].tick_params(which = 'both', direction='in')
axes[1,1].text(0.1*n4, 0.85*n4, r'$a/\lambda = 4$', backgroundcolor='white')


fig.colorbar(im4, ax = axes, location='bottom', shrink=0.97, aspect=40)



plt.savefig('PEC_Es_2D.svg'.format(c), format='svg')
plt.close()
