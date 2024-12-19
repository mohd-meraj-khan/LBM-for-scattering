import numpy as np
import sys
import os
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt

from scipy.signal import hilbert
from peakdetect import peakdetect


directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)


directory1 = 'gallery'
if not os.path.exists(directory1):
    os.makedirs(directory1)






ratio = 4

er2 = 2


r = [2, 3, 4, 5]



for j in range(len(r)):
    
    Escat = np.loadtxt(directory+'/Ez_scat{}_{}.txt'.format(r[j], ratio))

    RCS_LBM = []

    for i in range(len(Escat[0, :])):
        peaks = peakdetect(Escat[:, i], lookahead=10)
        higherPeaks = np.array(peaks[0])
        np.average(higherPeaks[:,1])
        RCS_LBM.append(2 * np.pi * r[j] * ratio * np.average(higherPeaks[:,1])**2)

    ###################################

    brcs = open(directory+"/BRCS_LBM{}_{}.txt".format(r[j], ratio), "w")
    np.savetxt(brcs, RCS_LBM)
    brcs.close()

    MRCS = np.array([RCS_LBM[180]])

    mrcs = open(directory+"/MRCS_LBM{}_{}.txt".format(r[j], ratio), "w")
    np.savetxt(mrcs, MRCS)
    mrcs.close()

    ###################################



    ###################################
    plt.clf()

    plt.rc('font', family = 'serif', size = 10)
    plt.rc('xtick', labelsize = 10)
    plt.rc('ytick', labelsize = 10)
    plt.rc('lines', markersize = 2, lw = 0.75)
    plt.rc('text', usetex = True)
    ######################################

    n = 100
    
    ################################################################################
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 1.75), dpi=600, constrained_layout = True)

    peaks = peakdetect(Escat[:, n], lookahead=10)
    higherPeaks = np.array(peaks[0])
    lowerPeaks = np.array(peaks[1])
    
    print(np.absolute(np.average(higherPeaks[:,1])))
    print(np.absolute(np.average(lowerPeaks[:,1])))
    
    plt.plot(np.arange(len(Escat[:, n])), Escat[:, n] , 'k-')
    plt.plot(higherPeaks[:,0], higherPeaks[:,1], 'ro')
    plt.plot(lowerPeaks[:,0], lowerPeaks[:,1], 'ko')

    plt.grid()
        
    plt.savefig(directory1+'/EzS_ratio_{}_r_{}.svg'.format(ratio, r[j]))
    plt.close()
    ################################################################################

