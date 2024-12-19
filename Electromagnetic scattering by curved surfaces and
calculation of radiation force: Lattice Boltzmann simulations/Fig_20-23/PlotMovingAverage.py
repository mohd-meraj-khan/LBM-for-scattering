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

er = 4.0
a = 100





ratio91 = 0.91
period91 = int(np.round(3 * a / ratio91))

ratio93 = 0.93
period93 = int(np.round(3 * a / ratio93))

ratio95 = 0.95
period95 = int(np.round(3 * a / ratio95))


U91  = np.loadtxt('energy_0.91.txt')
U93  = np.loadtxt('energy_0.93.txt')
U95  = np.loadtxt('energy_0.95.txt')

fx91 = np.loadtxt('fxSum_ratio_0.91.txt')
fx93 = np.loadtxt('fxSum_ratio_0.93.txt')
fx95 = np.loadtxt('fxSum_ratio_0.95.txt')


maForce91 = []
maEnergy91 = []

maForce93 = []
maEnergy93 = []

maForce95 = []
maEnergy95 = []



step = 1

U = 10000


i = 0

while (i < len(fx91) - period91 + 1):
    window = fx91[i : i + period91]
    window_average = sum(window) / period91
    maForce91.append(window_average)
    i += step

i = 0

while (i < len(U91) - period91 + 1):
    window = U91[i : i + period91]
    window_average = sum(window) / period91
    maEnergy91.append(window_average / U)
    i += step

xmaForce91 = np.arange(0, len(maForce91))/period91
xmaEnergy91 = np.arange(0, len(maEnergy91))/period91

#################################

i = 0

while (i < len(fx93) - period93 + 1):
    window = fx93[i : i + period93]
    window_average = sum(window) / period93
    maForce93.append(window_average)
    i += step

i = 0

while (i < len(U93) - period93 + 1):
    window = U93[i : i + period93]
    window_average = sum(window) / period93
    maEnergy93.append(window_average / U)
    i += step

xmaForce93 = np.arange(0, len(maForce93))/period93
xmaEnergy93 = np.arange(0, len(maEnergy93))/period93


###################################

i = 0

while (i < len(fx95) - period95 + 1):
    window = fx95[i : i + period95]
    window_average = sum(window) / period95
    maForce95.append(window_average)
    i += step

i = 0

while (i < len(U95) - period95 + 1):
    window = U95[i : i + period95]
    window_average = sum(window) / period95
    maEnergy95.append(window_average / U)
    i += step

xmaForce95 = np.arange(0, len(maForce95))/period95
xmaEnergy95 = np.arange(0, len(maEnergy95))/period95
###############################


N = 239000
xma91 = np.arange(len(maForce91[:N]), len(maForce91), 1)/period91
xma93 = np.arange(len(maForce93[:N]), len(maForce93), 1)/period93
xma95 = np.arange(len(maForce95[:N]), len(maForce95), 1)/period95

##M = 220000
x91 = np.arange(len(fx91[:N]), len(fx91), 1)/period91
x93 = np.arange(len(fx93[:N]), len(fx93), 1)/period93
x95 = np.arange(len(fx95[:N]), len(fx95), 1)/period95

print(len(maForce91[:N]), len(maForce91))
print(len(maForce93[:N]), len(maForce93))
print(len(maForce95[:N]), len(maForce95))

print(len(fx91[:N]), len(fx91))
print(len(fx93[:N]), len(fx93))
print(len(fx95[:N]), len(fx95))

    


##############################################################################################################################################################################


plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)




##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2.15), dpi=600, constrained_layout = True)


ax.plot(xmaEnergy91, maEnergy91, 'k-')
ax.plot(xmaEnergy93, maEnergy93, 'r-')
ax.plot(xmaEnergy95, maEnergy95, 'b-')


ax.set_xlabel(r'$t/T$')
ax.set_ylabel(r'$\frac{\left< U \right> / L}{\varepsilon_0 E_0^2 a^2}$')

plt.legend([r'$a/ \lambda = 0.91$', r'$a/ \lambda = 0.93$', r'$a/ \lambda = 0.95$'])


##plt.grid()

plt.savefig('energy_er_{}.svg'.format(er))
plt.close(fig)

##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2.15), dpi=600, constrained_layout = True)


ax.plot(xmaForce91, maForce91, 'k-')
plt.axhline(y=0.521852, color='k', linestyle='--')
ax.plot(xmaForce93, maForce93, 'r-')
plt.axhline(y=0.620928, color='r', linestyle='--')
ax.plot(xmaForce95, maForce95, 'b-')
plt.axhline(y=1.44645, color='b', linestyle='--')

ax.set_xlabel(r'$t/T$')
ax.set_ylabel(r'$\frac{\left< F_x \right> / L}{\lambda \varepsilon_0 E_0^2}$')


R91 = patches.Patch(color='k', label=r'$a/ \lambda = 0.91$')
R93 = patches.Patch(color='r', label=r'$a/ \lambda = 0.93$')
R95 = patches.Patch(color='b', label=r'$a/ \lambda = 0.95$')



legend = plt.legend(ncol=1, handles=[R91, R93, R95], loc='best', bbox_to_anchor=(0.5, 0.4, 0.5, 0.5))
ax.add_artist(legend)


##plt.grid()

plt.savefig('MA_er_{}.svg'.format(er))
plt.close(fig)
##################################################################################################


##################################################################################################

fig, ax = plt.subplots(figsize = (2.69, 1.85), dpi=600, constrained_layout = True)


plt.axhline(y=0.521852, color='k', linestyle='--')
ax.plot(x91, fx91[N:], 'b-')
ax.plot(xma91, maForce91[N:], 'r-')

ax.set_xlabel(r'$t/T$')
ax.set_ylabel(r'$\frac{ F_x  / L}{\lambda \varepsilon_0 E_0^2}$')

plt.legend([r'$ \langle $ Analytical $ \rangle $', r'LBM', r'$ \langle $ LBM $ \rangle $'])


plt.savefig('MA_er_{}_0.91.svg'.format(er))
plt.close(fig)
##################################################################################################


##################################################################################################

fig, ax = plt.subplots(figsize = (2.69, 1.85), dpi=600, constrained_layout = True)


plt.axhline(y=0.620928, color='k', linestyle='--')
ax.plot(x93, fx93[N:], 'b-')
ax.plot(xma93, maForce93[N:], 'r-')

ax.set_xlabel(r'$t/T$')
##ax.set_ylabel(r'$\frac{ F_x  / L}{\lambda \varepsilon_0 E_0^2}$')


plt.savefig('MA_er_{}_0.93.svg'.format(er))
plt.close(fig)
##################################################################################################


##################################################################################################

fig, ax = plt.subplots(figsize = (2.69, 1.85), dpi=600, constrained_layout = True)


plt.axhline(y=1.44645, color='k', linestyle='--')
ax.plot(x95, fx95[N:], 'b-')
ax.plot(xma95, maForce95[N:], 'r-')

ax.set_xlabel(r'$t/T$')
##ax.set_ylabel(r'$\frac{ F_x  / L}{\lambda \varepsilon_0 E_0^2}$')


plt.savefig('MA_er_{}_0.95.svg'.format(er))
plt.close(fig)
##################################################################################################



##################################################################################################

fig = plt.figure(figsize = (6.69, 1.85), dpi=600)

gs=GridSpec(1,3)


ax0 = fig.add_subplot(gs[0,0])
plt.axhline(y=0.521852, color='k', linestyle='--')
ax0.plot(x91, fx91[N:], 'b-')
ax0.plot(xma91, maForce91[N:], 'r-')

plt.text(724.25, 1.42, r'(a)', backgroundcolor='white')

##plt.legend([r'$ \langle $ Analytical $ \rangle $', r'LBM', r'$ \langle $ LBM $ \rangle $'])

ax0.set_xlabel(r'$t/T$')
ax0.set_ylabel(r'$\frac{ F_x  / L}{\lambda \varepsilon_0 E_0^2}$')




ax1 = fig.add_subplot(gs[0,1])
plt.axhline(y=0.620928, color='k', linestyle='--')
ax1.plot(x93, fx93[N:], 'b-')
ax1.plot(xma93, maForce93[N:], 'r-')

plt.text(740.25, 1.32, r'(b)', backgroundcolor='white')

ax1.set_xlabel(r'$t/T$')
##ax.set_ylabel(r'$\frac{ F_x  / L}{\lambda \varepsilon_0 E_0^2}$')



ax2 = fig.add_subplot(gs[0,2])
plt.axhline(y=1.44645, color='k', linestyle='--')
ax2.plot(x95, fx95[N:], 'b-')
ax2.plot(xma95, maForce95[N:], 'r-')

plt.text(756.5, 1.85, r'(c)', backgroundcolor='white')

ax2.set_xlabel(r'$t/T$')
##ax2.set_ylabel(r'$\frac{ F_x  / L}{\lambda \varepsilon_0 E_0^2}$')

plt.tight_layout()

plt.savefig('F_ins_Avg_er_{}.svg'.format(er))
plt.close(fig)
##################################################################################################


