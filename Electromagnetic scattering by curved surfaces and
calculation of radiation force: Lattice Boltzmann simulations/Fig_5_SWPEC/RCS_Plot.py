import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import sys
import os


directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)


ratio = [0.05, 0.1, 0.5, 1, 2, 4]

dphi = 1
phi = np.arange(0, 180 + 0.001, dphi)

dR = 0.1
R = np.arange(2, 5 + 0.001, dR)
RLBM = [2, 3, 4, 5]



BRCS_LBM = []
MRCS_LBM = np.zeros((len(ratio), len(RLBM)))

BRCS_exact = []
MRCS_exact = []


for i in range (len(ratio)):
    for j in range(len(RLBM)):
        MRCS_LBM[i, j] = np.loadtxt(directory+'/MRCS_LBM{}_{}.txt'.format(RLBM[j], ratio[i]))


for i in range (len(ratio)):
    BRCS_LBM.append(np.loadtxt(directory+'/BRCS_LBM5_{}.txt'.format(ratio[i])))
    
    BRCS_exact.append(np.loadtxt(directory+'/BRCS_exact_{}.txt'.format(ratio[i])))
    MRCS_exact.append(np.loadtxt(directory+'/MRCS_exact_{}.txt'.format(ratio[i])))




###################################


plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)

################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 4.75), dpi=600, constrained_layout = True)

ax1 = ax.inset_axes([0.25, 0.7, 0.71, 0.27])

############ BRCS ############
plt.plot(phi, BRCS_exact[0], 'k-')
plt.plot(phi, BRCS_exact[1], 'r-')
plt.plot(phi, BRCS_exact[2], 'b-')
plt.plot(phi, BRCS_exact[3], 'c-')
plt.plot(phi, BRCS_exact[4], 'm-')
plt.plot(phi, BRCS_exact[5], color='dimgrey', linestyle='-')

plt.plot(phi, BRCS_LBM[0], 'k--')
plt.plot(phi, BRCS_LBM[1], 'r--')
plt.plot(phi, BRCS_LBM[2], 'b--')
plt.plot(phi, BRCS_LBM[3], 'c--')
plt.plot(phi, BRCS_LBM[4], 'm--')
plt.plot(phi, BRCS_LBM[5], color='dimgrey', linestyle='--')

######### MRCS ############
ax1.plot(R, MRCS_exact[0], 'k-')
ax1.plot(R, MRCS_exact[1], 'r-')
ax1.plot(R, MRCS_exact[2], 'b-')
ax1.plot(R, MRCS_exact[3], 'c-')
ax1.plot(R, MRCS_exact[4], 'm-')
ax1.plot(R, MRCS_exact[5], color='dimgrey', linestyle='-')

ax1.plot(RLBM, MRCS_LBM[0], 'ko')
ax1.plot(RLBM, MRCS_LBM[1], 'ro')
ax1.plot(RLBM, MRCS_LBM[2], 'bo')
ax1.plot(RLBM, MRCS_LBM[3], 'co')
ax1.plot(RLBM, MRCS_LBM[4], 'mo')
ax1.plot(RLBM, MRCS_LBM[5], color='dimgrey',  linestyle='none', marker='o')

ax1.set_xlabel(r'$r/a$')

ax1.set_yscale('log')
###########################

plt.xticks([0, 30, 60, 90, 120, 150, 180])

ra005 = patches.Patch(color='k', label=r'$a/ \lambda = 0.05$')
ra01  = patches.Patch(color='r', label=r'$a/ \lambda = 0.1$')
ra05  = patches.Patch(color='b', label=r'$a/ \lambda = 0.5$')
ra1   = patches.Patch(color='c', label=r'$a/ \lambda = 1$')
ra2   = patches.Patch(color='m', label=r'$a/ \lambda = 2$')
ra4   = patches.Patch(color='dimgrey', label=r'$a/ \lambda = 4$')

legend = ax.legend(ncol=2, mode='expand', loc=[0.02, 1.02], handles=[ra005, ra01, ra05, ra1, ra2, ra4])#, , bbox_to_anchor=(0.4, 1.2), )


plt.xlabel(r'$\phi (^o)$')
plt.ylabel(r'$\sigma / \lambda$')


plt.yscale('log')

    
plt.savefig('RCS_PEC.svg')
plt.close()
################################################################################





################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 1.75), dpi=600, constrained_layout = True)


plt.plot(phi, np.absolute((BRCS_exact[0] - BRCS_LBM[0]) / BRCS_exact[0] * 100) , 'k-')
plt.plot(phi, np.absolute((BRCS_exact[1] - BRCS_LBM[1]) / BRCS_exact[1] * 100), 'r-')
plt.plot(phi, np.absolute((BRCS_exact[2] - BRCS_LBM[2]) / BRCS_exact[2] * 100), 'b-')
plt.plot(phi, np.absolute((BRCS_exact[3] - BRCS_LBM[3]) / BRCS_exact[3] * 100), 'c-')
plt.plot(phi, np.absolute((BRCS_exact[4] - BRCS_LBM[4]) / BRCS_exact[4] * 100), 'm-')
plt.plot(phi, np.absolute((BRCS_exact[5] - BRCS_LBM[5]) / BRCS_exact[5] * 100), color='dimgrey', linestyle='-')


plt.xticks([0, 30, 60, 90, 120, 150, 180])


plt.xlabel(r'$\phi (^o)$')
plt.ylabel(r'$\%$ Error')

    
plt.savefig('RCS_PEC_Err.svg')
plt.close()
################################################################################

