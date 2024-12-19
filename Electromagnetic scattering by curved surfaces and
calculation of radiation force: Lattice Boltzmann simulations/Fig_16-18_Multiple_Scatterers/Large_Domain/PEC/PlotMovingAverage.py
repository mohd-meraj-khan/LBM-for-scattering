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

er = 10000
a = 40





ratio2 = 2
period2 = int(np.round(3 * a / ratio2))


U2  = np.loadtxt('data/energy_ratio_1_a_40.txt')



maEnergy2 = []





step = 1

U = 1600




i = 0

while (i < len(U2) - period2 + 1):
    window = U2[i : i + period2]
    window_average = sum(window) / period2
    maEnergy2.append(window_average / U)
    i += step

xmaEnergy2 = np.arange(0, len(maEnergy2))/period2

#################################




###################################


###############################


##N = 239000
##xma91 = np.arange(len(maForce91[:N]), len(maForce91), 1)/period91
##xma93 = np.arange(len(maForce93[:N]), len(maForce93), 1)/period93
##xma95 = np.arange(len(maForce95[:N]), len(maForce95), 1)/period95
##
####M = 220000
##x91 = np.arange(len(fx91[:N]), len(fx91), 1)/period91
##x93 = np.arange(len(fx93[:N]), len(fx93), 1)/period93
##x95 = np.arange(len(fx95[:N]), len(fx95), 1)/period95
##
##print(len(maForce91[:N]), len(maForce91))
##print(len(maForce93[:N]), len(maForce93))
##print(len(maForce95[:N]), len(maForce95))
##
##print(len(fx91[:N]), len(fx91))
##print(len(fx93[:N]), len(fx93))
##print(len(fx95[:N]), len(fx95))

    


##############################################################################################################################################################################


plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)




##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2.15), dpi=600, constrained_layout = True)


ax.plot(xmaEnergy2, maEnergy2, 'k-')
##ax.plot(xmaEnergy93, maEnergy93, 'r-')
##ax.plot(xmaEnergy95, maEnergy95, 'b-')


ax.set_xlabel(r'$t/T$')
ax.set_ylabel(r'$\frac{\left< U \right> / L}{\varepsilon_0 E_0^2 a^2}$')

plt.legend([r'$a/ \lambda = 0.91$', r'$a/ \lambda = 0.93$', r'$a/ \lambda = 0.95$'])


##plt.grid()

plt.savefig('energy_er_{}.svg'.format(er))
plt.close(fig)

##################################################################################################




