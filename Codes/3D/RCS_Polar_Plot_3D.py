import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import sys
import os

from Module_Parameters_3D import *


directory = 'data/rcs'



plots = 'plots/plots_RCS/polar'
if not os.path.exists(plots):
    os.makedirs(plots)




phi = np.linspace(0, 2*np.pi, 361)




theta = np.linspace(0, 2*np.pi, 361)



BRCS_LBM   = np.load(directory+'/RCS_LBM_{}_{}_{}.npy'.format(ratio, er2, phi0))
BRCS_exact = np.load(directory+'/BRCS_exact_{}_{}_{}.npy'.format(ratio, er2, phi0))

    




###################################
width = 2.23

plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)

################################################################################
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize = (width, width), dpi=600, constrained_layout = True)




ax.plot(theta, BRCS_exact, 'k-')

ax.plot(theta, BRCS_LBM, 'r--')

ax.set_rscale('symlog')

ax.set_rlabel_position(-125)
ax.grid(True)

ax.set_thetagrids([45, 135, 225, 315], labels=[r'$\pi/4$', r'$3\pi/4$', r'$5\pi/4$', r'$7\pi/4$'])


##plt.legend([r'Analytical', r'LBM'])


    
plt.savefig(plots+'/RCS_Polar_{}_{}_{}.svg'.format(ratio, er2, phi0))
plt.close()
################################################################################




################################################################################
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize = (width, width), dpi=600, constrained_layout = True)




ax.plot(theta, np.absolute((BRCS_exact - BRCS_LBM) / BRCS_exact)*100, 'k-')



ax.set_rscale('symlog')

ax.set_rlabel_position(-125)
ax.grid(True)

ax.set_thetagrids([45, 135, 225, 315], labels=[r'$\pi/4$', r'$3\pi/4$', r'$5\pi/4$', r'$7\pi/4$'])


    
plt.savefig(plots+'/RCS_Polar_{}_{}_{}_Err.svg'.format(ratio, er2, phi0))
plt.close()
################################################################################

print(np.shape(BRCS_LBM))
print(np.shape(BRCS_exact))
print(np.shape(np.absolute((BRCS_exact - BRCS_LBM) / BRCS_exact)*100))



