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


##A = 0.5
##
##a1  = 50
##Nx1 = 200
##Ny1 = 200
##
##a2  = 100
##Nx2 = 400
##Ny2 = 400
##
##a3  = 150
##Nx3 = 600
##Ny3 = 600


Mitri05 = np.loadtxt('Mitri_0.5_corrugated.txt')
Mitri10 = np.loadtxt('Mitri_1.0_corrugated.txt')
Mitri15 = np.loadtxt('Mitri_1.5_corrugated.txt')

fx05 = np.loadtxt('Fx_A_0.5.txt')
fx10 = np.loadtxt('Fx_A_1.txt')
fx15 = np.loadtxt('Fx_A_1.5.txt')


##############################################################################################################################################################################



x_new = np.linspace(0.04, 0.76, 13)*2*np.pi


x_Mitri05 = Mitri05[:, 0]
y_Mitri05 = Mitri05[:, 1]

x_Mitri10 = Mitri10[:, 0]
y_Mitri10 = Mitri10[:, 1]

x_Mitri15 = Mitri15[:, 0]
y_Mitri15 = Mitri15[:, 1]

f05 = interpolate.interp1d(x_Mitri05, y_Mitri05)
f10 = interpolate.interp1d(x_Mitri10, y_Mitri10)
f15 = interpolate.interp1d(x_Mitri15, y_Mitri15)

y_Mitri05_new = f05(x_new)
y_Mitri10_new = f10(x_new)
y_Mitri15_new = f15(x_new)


err05 = np.absolute((y_Mitri05_new - fx05[:])) * 100
err10 = np.absolute((y_Mitri10_new - fx10[:])) * 100
err15 = np.absolute((y_Mitri15_new - fx15[:])) * 100


############################################################################################################################################################################

ratio = np.linspace(0.04, 0.76, 13)*2*np.pi


plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)



##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2.45), dpi=600, constrained_layout = True)




ax.plot(Mitri05[:,0], Mitri05[:,1], 'k-', ratio, fx05, 'ko')
ax.plot(Mitri10[:,0], Mitri10[:,1], 'r-', ratio, fx10, 'ro')
ax.plot(Mitri15[:,0], Mitri15[:,1], 'b-', ratio, fx15, 'bo')


ax.legend([r'Mitri 2019', r'LBM ($a = 50 \Delta x$)', r'LBM ($a = 100 \Delta x$)', r'LBM ($a = 150 \Delta x$)'])
ax.set_xlabel(r'$k b$')
ax.set_ylabel(r'$\frac{\left< F_x \right>/L}{b \varepsilon_0 E_0^2}$')

a05  = patches.Patch(color='k', label=r'$A = 0.5$')
a10  = patches.Patch(color='r', label=r'$A = 1.0$')
a15  = patches.Patch(color='b', label=r'$A = 1.5$')


legend = plt.legend(handles=[a05, a10, a15])
ax.add_artist(legend)

##plt.text(2.5, 1.2, r'$L=3 a$')


##plt.grid()

plt.savefig('corrugatedPEC.svg')
plt.close(fig)


##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2.15), dpi=600, constrained_layout = True)


ax.plot(x_new, err05, 'k-o', x_new, err10, 'r-o', x_new, err15, 'b-o')
ax.legend([r'LBM ($a = 50 \Delta x$)', r'LBM ($a = 100 \Delta x$)', r'LBM ($a = 150 \Delta x$)'])
ax.set_xlabel(r'$k b$')
ax.set_ylabel(r'$\%$ Error')

a05  = patches.Patch(color='k', label=r'$A = 0.5$')
a10  = patches.Patch(color='r', label=r'$A = 1.0$')
a15  = patches.Patch(color='b', label=r'$A = 1.5$')


legend = plt.legend(handles=[a05, a10, a15])
ax.add_artist(legend)

##plt.text(3, 6, r'$L=3 a$')


plt.grid()

plt.savefig('Mitri_Fx_err.svg')
plt.close(fig)

