import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import sys
import os


directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)


er2 = 4

dratio = 0.001
ratio = np.arange(0.9, 1.1 + dratio, dratio)

ratio_LBM = np.array([0.90, 0.925, 0.95, 0.975, 1.0, 1.025, 1.05, 1.075, 1.10])




Q_exact = np.loadtxt(directory+'/Q_sca_exact_{}.txt'.format(er2))
Q_LBM   = np.loadtxt(directory+'/Q_sca_{}.txt'.format(er2))





###################################


plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)

################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.15), dpi=600, constrained_layout = True)



############ BRCS ############
plt.plot(ratio, Q_exact, 'k-')
plt.plot(ratio_LBM, Q_LBM, 'ro')


##plt.xticks([0, 30, 60, 90, 120, 150, 180])



plt.xlabel(r'$a / \lambda$')
plt.ylabel(r'$C_{sca} / (2 a)$')

plt.legend([r'Analytical', r'LBM']) 

##plt.grid()

##plt.xscale('log')

    
plt.savefig('Q_sca_{}_Mie.svg'.format(er2))
plt.close()
################################################################################


