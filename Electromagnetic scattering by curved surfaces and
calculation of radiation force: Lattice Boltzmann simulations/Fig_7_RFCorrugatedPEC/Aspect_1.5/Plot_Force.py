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


path = 'data'

A = 1.5

Mitri = np.loadtxt('Mitri_1.5_corrugated.txt')

FxRayleigh = np.loadtxt('FxAvgRayleigh_A_{}.txt'.format(A))
FxMie = np.loadtxt('FxAvgMie_A_{}.txt'.format(A))

Fx = np.concatenate((FxRayleigh, FxMie))



fx = open("Fx_A_{}.txt".format(A), "w")
np.savetxt(fx, Fx)
fx.close()

##############################################################################################################################################################################


x_new = np.linspace(0.04, 0.76, 13)*2*np.pi


x_Mitri = Mitri[:, 0]
y_Mitri = Mitri[:, 1]

f = interpolate.interp1d(x_Mitri, y_Mitri)

y_Mitri_new = f(x_new)


err = np.absolute((y_Mitri_new - Fx)) * 100






############################################################################################################################################################################

ratio = np.linspace(0.04, 0.76, 13)*2*np.pi


plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 1)
plt.rc('text', usetex = True)



##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2.15), dpi=600, constrained_layout = True)




ax.plot(Mitri[:,0], Mitri[:,1], 'k-', ratio, Fx, 'ro')
ax.legend([r'Mitri 2019', r'LBM'])
ax.set_xlabel(r'$2 \pi b / \lambda$')
ax.set_ylabel(r'$\frac{F_x/L}{b \varepsilon_0 E_0^2}$')

##plt.text(2.5, 1.2, r'$L=3 a$')


##plt.grid()

plt.savefig('Mitri_Fx_aspect_1.0.svg')
plt.close(fig)


##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2.15), dpi=600, constrained_layout = True)


ax.plot(x_new, err, 'r-o')
##ax.legend([r'LBM ($a = 50 \Delta x$)'])
ax.set_xlabel(r'$2 \pi b / \lambda$')
ax.set_ylabel(r'$\%$ Error')

##plt.text(3, 6, r'$L=3 a$')


plt.grid()

plt.savefig('Mitri_Fx_aspect_1.0_err.svg')
plt.close(fig)

