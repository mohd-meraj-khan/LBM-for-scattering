import numpy as np
import scipy.special as sc

import sys
import os
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt

from scipy.signal import hilbert, find_peaks
from peakdetect import peakdetect



ratio = np.array([0.90, 0.925, 0.95, 0.975, 1.0, 1.025, 1.05, 1.075, 1.10])
a = 100

n = 4
Nx, Ny, Nz = n*a, n*a, 1

er1, mur1, er2 = 1, 1, 4


dtheta = 1
theta = np.arange(0, 360 + 0.01, dtheta)

# surface normal (the integration is done at a circle enclosing the scatterer)
nx = np.cos(theta*np.pi/180)
ny = np.sin(theta*np.pi/180)
nz = 0


R = 1.0*a


directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)


directory1 = 'gallery'
if not os.path.exists(directory1):
    os.makedirs(directory1)





velocity = 1.0/3

E0 = 1

I_inc = 0.5 * np.sqrt(er1 / mur1) * E0**2




Q_sca = np.zeros(len(ratio))



for m in range(len(ratio)):

    ######### DON'T CHANGE ########
    wavelength = a / ratio[m]
    period = wavelength / velocity
    omega = 2 * np.pi / period
    ###############################

    k = 2 * np.pi / wavelength
        
    # loading the field data
    Ez_scat = np.loadtxt(directory+'/Ez_scat_{}.txt'.format(ratio[m]))
    Hx_scat = np.loadtxt(directory+'/Hx_scat_{}.txt'.format(ratio[m]))
    Hy_scat = np.loadtxt(directory+'/Hy_scat_{}.txt'.format(ratio[m]))
    


    # initializing the frequency domain fields
    Hx_phasor = np.zeros(len(Ez_scat[0, :]), dtype=complex)
    Hy_phasor = np.zeros(len(Ez_scat[0, :]), dtype=complex)
    Ez_phasor = np.zeros(len(Ez_scat[0, :]), dtype=complex)



    # computing the frequency domain currents
    Sum_Hx_phasor = 0
    Sum_Hy_phasor = 0
    Sum_Ez_phasor = 0




    for i in range(len(Ez_scat[:, 0])):

        Sum_Hx_phasor += Hx_scat[i, :] * np.exp((0 - 1j) * omega * i) / len(Hx_scat[:, 0]) * 2
        Sum_Hy_phasor += Hy_scat[i, :] * np.exp((0 - 1j) * omega * i) / len(Hy_scat[:, 0]) * 2
        Sum_Ez_phasor += Ez_scat[i, :] * np.exp((0 - 1j) * omega * i) / len(Ez_scat[:, 0]) * 2

       
    Hx_phasor = Sum_Hx_phasor
    Hy_phasor = Sum_Hy_phasor
    Ez_phasor = Sum_Ez_phasor



    Hr_phasor   =   Hx_phasor * nx + Hy_phasor * ny
    Hphi_phasor = - Hx_phasor * ny + Hy_phasor * nx


    W_sca = - 0.5 * np.sum(np.real(Ez_phasor * np.conjugate(Hphi_phasor))) * R * (dtheta * np.pi/180)

    C_sca = W_sca / I_inc

    Q_sca[m] = C_sca / (2*a)

    print(Q_sca)


##    RCS_LBM_far.append(np.absolute(rcs))

##print(RCS_LBM_far)


# calculation of far-field bistatic SW
Qsca = open(directory+"/Q_sca_{}.txt".format(er2), "w")
np.savetxt(Qsca, Q_sca)
Qsca.close()
###################################

sys.exit()

n = 100


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
plt.plot(np.arange(len(Ez_scat[:, n])), Ez_scat[:, n], 'k-')
plt.plot(np.arange(len(Hx_scat[:, n])), Hx_scat[:, n], 'r-')
plt.plot(np.arange(len(Hy_scat[:, n])), Hy_scat[:, n], 'b-')


##plt.xticks([0, 30, 60, 90, 120, 150, 180])

plt.legend([r'$E_z$', r'$H_x$', r'$H_y$'])

plt.ylabel(r'Field')
plt.xlabel(r'time')

##plt.ylim(2, 8)

plt.grid()

##plt.xscale('log')

    
plt.savefig('Field_{}.svg'.format(er2))
plt.close()
################################################################################


