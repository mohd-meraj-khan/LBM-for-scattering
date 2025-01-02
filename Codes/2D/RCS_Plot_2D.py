import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import sys
import os

from Module_Parameters_2D import *


directory = 'data/rcs'



plots = 'plots/plots_RCS'
if not os.path.exists(plots):
    os.makedirs(plots)




dphi = 1
phi = np.arange(0, 360 + 0.001, dphi)



BRCS_LBM = np.load(directory+'/RCS_LBM_{}.npy'.format(ratio))
BRCS_exact = np.load(directory+'/BRCS_exact_{}.npy'.format(ratio))
    




###################################


plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)

################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.05), dpi=600, constrained_layout = True)




plt.plot(phi, BRCS_exact, 'k-')

plt.plot(phi, BRCS_LBM, 'r--')


plt.xticks([0, 60, 120, 180, 240, 300, 360])

plt.xlabel(r'$\phi (^o)$')
plt.ylabel(r'$\sigma / \lambda$')

plt.grid(which='both')


plt.yscale('log')

    
plt.savefig(plots+'/RCS_{}.svg'.format(ratio))
plt.close()
################################################################################







