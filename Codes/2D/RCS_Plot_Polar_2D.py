import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import sys
import os

from Module_Parameters_2D import *


directory = 'data/rcs'



plots = 'plots/plots_RCS/polar'
if not os.path.exists(plots):
    os.makedirs(plots)



phi = np.linspace(0, 2*np.pi, 361)


BRCS_LBM = np.load(directory+'/RCS_LBM_{}_{}.npy'.format(er2, ratio))
BRCS_exact = np.load(directory+'/BRCS_exact_{}_{}.npy'.format(er2, ratio))
    




###################################


plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)

################################################################################
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize = (3.35*0.859, 3.35*0.859), dpi=600, constrained_layout = True)



ax.plot(phi, BRCS_exact, 'k-')

ax.plot(phi, BRCS_LBM, 'r--')


ax.set_rscale('symlog')

ax.set_rlabel_position(-125)
ax.grid(True)

    
plt.savefig(plots+'/RCS_{}_{}.svg'.format(er2, ratio))
plt.close()
################################################################################




################################################################################
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize = (3.35*0.859, 3.35*0.859), dpi=600, constrained_layout = True)


ax.plot(phi, np.absolute((BRCS_exact - BRCS_LBM) / BRCS_exact)*100, 'k.')


ax.set_rlabel_position(-125)
ax.grid(True)

    
plt.savefig(plots+'/RCS_{}_{}_Err.svg'.format(er2, ratio))
plt.close()
################################################################################



