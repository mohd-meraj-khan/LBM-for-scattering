import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

import sys
import os


directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)


er2 = [2, 5, 10, 20]

ratio = 0.5

dphi = 1
phi = np.arange(0, 180 + 0.001, dphi)

dR = 0.1
R = np.arange(2, 5 + 0.001, dR)
RLBM = [2, 3, 4, 5]



BRCS_LBM = []
MRCS_LBM = np.zeros((len(er2), len(RLBM)))

BRCS_exact = []
MRCS_exact = []


for i in range (len(er2)):
    for j in range(len(RLBM)):
        MRCS_LBM[i, j] = np.loadtxt(directory+'/MRCS_LBM{}_{}.txt'.format(RLBM[j], er2[i]))


for i in range (len(er2)):
    BRCS_LBM.append(np.loadtxt(directory+'/BRCS_LBM5_{}.txt'.format(er2[i])))
    
    BRCS_exact.append(np.loadtxt(directory+'/BRCS_exact_{}.txt'.format(er2[i])))
    MRCS_exact.append(np.loadtxt(directory+'/MRCS_exact_{}.txt'.format(er2[i])))

##MRCS_exact.append(np.loadtxt(directory+'/MRCS_exact_PEC.txt'))


###################################


plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)

################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2), dpi=600, constrained_layout = True)


############ BRCS ############
plt.plot(phi, BRCS_exact[0], 'b-')
plt.plot(phi, BRCS_exact[1], 'r-')
plt.plot(phi, BRCS_exact[2], 'c-')
plt.plot(phi, BRCS_exact[3], 'k-')

plt.plot(phi, BRCS_LBM[0], 'b--')
plt.plot(phi, BRCS_LBM[1], 'r--')
plt.plot(phi, BRCS_LBM[2], 'c--')
plt.plot(phi, BRCS_LBM[3], 'k--')




plt.xticks([0, 30, 60, 90, 120, 150, 180])

ra2 = patches.Patch(color='b', label=r'$\varepsilon_r = 2$')
ra5  = patches.Patch(color='r', label=r'$\varepsilon_r = 5$')
ra10  = patches.Patch(color='c', label=r'$\varepsilon_r = 10$')
ra20  = patches.Patch(color='k', label=r'$\varepsilon_r = 20$')
##rapec  = patches.Patch(color='k', label=r'PEC')


legend = ax.legend(ncol=2, handles=[ra2, ra5, ra10, ra20], loc='upper right')


plt.xlabel(r'$\phi (^o)$')
plt.ylabel(r'$\sigma / \lambda$')

plt.ylim(0.005, 300)
plt.yscale('log')

    
plt.savefig('BRCS_ratio_0.5.svg')
plt.close()
################################################################################



################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2), dpi=600, constrained_layout = True)



######### MRCS ############
ax.plot(R, MRCS_exact[0], 'b-')
ax.plot(R, MRCS_exact[1], 'r-')
ax.plot(R, MRCS_exact[2], 'c-')
ax.plot(R, MRCS_exact[3], 'k-')


ax.plot(RLBM, MRCS_LBM[0], 'bo')
ax.plot(RLBM, MRCS_LBM[1], 'ro')
ax.plot(RLBM, MRCS_LBM[2], 'co')
ax.plot(RLBM, MRCS_LBM[3], 'ko')

ax.set_xlabel(r'$r/a$')
plt.ylabel(r'$\sigma / \lambda$')

plt.yscale('log')

plt.savefig('MRCS_ratio_0.5.svg')
plt.close()
################################################################################





################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 1.75), dpi=600, constrained_layout = True)


plt.plot(phi, np.absolute((BRCS_exact[0] - BRCS_LBM[0]) / BRCS_exact[0] * 100), 'b-')
plt.plot(phi, np.absolute((BRCS_exact[1] - BRCS_LBM[1]) / BRCS_exact[1] * 100), 'r-')
plt.plot(phi, np.absolute((BRCS_exact[2] - BRCS_LBM[2]) / BRCS_exact[2] * 100), 'c-')
plt.plot(phi, np.absolute((BRCS_exact[3] - BRCS_LBM[3]) / BRCS_exact[3] * 100), 'k-')


plt.xticks([0, 30, 60, 90, 120, 150, 180])


plt.xlabel(r'$\phi (^o)$')
plt.ylabel(r'$\%$ Error')

    
plt.savefig('RCS_ratio_0.5_Err.svg')
plt.close()
################################################################################

