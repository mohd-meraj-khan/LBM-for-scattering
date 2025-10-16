import numpy as np
import matplotlib.pyplot as plt

import matplotlib.patches as patches
import sys
import os

from Module_Parameters_2D import *




directory = 'data/moving_average'


force = 'plots/MA_Force'
if not os.path.exists(force):
    os.makedirs(force)

energy = 'plots/MA_Energy'
if not os.path.exists(energy):
    os.makedirs(energy)



u  = np.load(directory+"/energy_{}_{}.npy".format(er2, ratio))


maEnergy = []


step = 1

U = a**2




i = 0

while (i < len(u) - period + 1):
    window_u = u[i : i + period]
    window_u_average = sum(window_u) / period
    maEnergy.append(window_u_average / U)
    
    i += step
###################################

xma = np.arange(0, len(maEnergy))/period
x = np.arange(0, len(u))/period




##############################################################################################################################################################################


plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)




n = 0
##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2.15), dpi=600, constrained_layout = True)

ax.plot(xma[n:], maEnergy[n:], 'k-')

ax.set_xlabel(r'$t / T$')
ax.set_ylabel(r'$\frac{\left< U \right> / L}{\varepsilon_0 E_0^2 a^2}$')

##plt.xscale('log')
plt.grid()

plt.savefig(energy+'/MA_energy_{}_{}.svg'.format(er2, ratio))
plt.close(fig)

##################################################################################################


sys.exit()


##################################################################################################


fx = np.load(directory+"/FxIns_{}_{}.npy".format(er2, ratio))
fy = np.load(directory+"/FyIns_{}_{}.npy".format(er2, ratio))
tz = np.load(directory+"/TzIns_{}_{}.npy".format(er2, ratio))

mafx = []
mafy = []
matz = []



while (i < len(u) - period + 1):
    window_fx = fx[i : i + period]
    window_fx_average = sum(window_fx) / period
    mafx.append(window_fx_average)
    
    window_fy = fy[i : i + period]
    window_fy_average = sum(window_fy) / period
    mafy.append(window_fy_average)

    window_tz = tz[i : i + period]
    window_tz_average = sum(window_tz) / period
    matz.append(window_tz_average)

    
    i += step
###################################


##################################################################################################

fig, ax = plt.subplots(figsize = (3.35, 2.15), dpi=600, constrained_layout = True)

ax.plot(xma[n:], mafx[n:], 'k-')
ax.plot(xma[n:], mafy[n:], 'r-')
ax.plot(xma[n:], matz[n:], 'b-')

ax.set_xlabel(r'$t / T$')

plt.legend([r'$\frac{\left< F_x \right> / L}{\lambda \varepsilon_0 E_0^2}$',
            r'$\frac{\left< F_y \right> / L}{\lambda \varepsilon_0 E_0^2}$',
            r'$\frac{\left< T_z \right> / L}{\lambda^2 \varepsilon_0 E_0^2}$'])

##plt.xscale('log')
plt.grid()

plt.savefig(force+'/MA_Force_{}_{}.svg'.format(er2, ratio))
plt.close(fig)

##################################################################################################
