import numpy as np

import sys
import os

from scipy.interpolate import RectBivariateSpline

import math
import scipy.special as sc
from math import e
import cmath
from scipy.optimize import curve_fit


directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)







ratio = 4

er2 = 2


r = 5

er1 = 1
mur1 = 1
mur2 = 1

v1 = 1 / (3 * np.sqrt(er1 * mur1))
v2 = 1 / (3 * np.sqrt(er2 * mur2))

alpha = np.sqrt((er1 * mur2) / (er2 * mur1))

k1a = 2 * np.pi * ratio
k2a = 2 * np.pi * v1 / v2 * ratio
k1r = 2 * np.pi * ratio * r


########################################## BRCS ###################################
dphi = 1
phi = np.arange(0, 180 + 0.001, dphi)

BRCS = np.zeros(len(phi))

n = 1000

for i in range(len(phi)):
   
    ezs_sum   = 0
        
    for l in range(0, n+1, 1):
            
        num  = (alpha * sc.jvp(l, k1a) * sc.jv(l, k2a) - sc.jv(l, k1a) * sc.jvp(l, k2a))
        den  = (sc.hankel2(l, k1a) * sc.jvp(l, k2a) - alpha * sc.h2vp(l, k1a) * sc.jv(l, k2a))
            
##        num = - sc.jv(l, k1a)
##        den = sc.hankel2(l, k1a)

        Cn = num / den

        if (l == 0):
            ezs = (0 - 1j)**l * (Cn * sc.hankel2(l, k1r)) * np.cos(l*phi[i]*np.pi/180)
        else :
            ezs = 2 * (0 - 1j)**l * (Cn * sc.hankel2(l, k1r)) * np.cos(l*phi[i]*np.pi/180)
            
        ezs_sum   += ezs

        if (np.absolute(ezs) <= 1e-10 and np.absolute(np.cos(l*phi[i]*np.pi/180)) >= 1e-10):
##            print(l)
            break       
  
    BRCS[i] = 2 * np.pi * r * ratio * np.absolute(ezs_sum)**2

###################################

rcs = open(directory+"/BRCS_exact_{}.txt".format(ratio), "w")
np.savetxt(rcs, BRCS)
rcs.close()


##rcs = open(directory+"/BRCS_exact_PEC.txt", "w")
##np.savetxt(rcs, BRCS)
##rcs.close()

###################################################################################





########################################## MRCS ###################################

dR = 0.1
R = np.arange(2, 5 + 0.001, dR)

MRCS = np.zeros(len(R))

for i in range(len(R)):
    k1a = 2 * np.pi * ratio
    k2a = 2 * np.pi * v1 / v2 * ratio
    k1r = 2 * np.pi * ratio * R[i]

    ezs_sum   = 0
        
    for l in range(0, n+1, 1):
            
        num  = (alpha * sc.jvp(l, k1a) * sc.jv(l, k2a) - sc.jv(l, k1a) * sc.jvp(l, k2a))
        den  = (sc.hankel2(l, k1a) * sc.jvp(l, k2a) - alpha * sc.h2vp(l, k1a) * sc.jv(l, k2a))
            
##        num = - sc.jv(l, k1a)
##        den = sc.hankel2(l, k1a)

        Cn = num / den

        if (l == 0):
            ezs = (0 - 1j)**l * (Cn * sc.hankel2(l, k1r)) * np.cos(l*np.pi)
        else :
            ezs = 2 * (0 - 1j)**l * (Cn * sc.hankel2(l, k1r)) * np.cos(l*np.pi)
                
        ezs_sum   += ezs

        if (np.absolute(ezs) <= 1e-10):
##            print(l)
            break

    MRCS[i] = 2 * np.pi * R[i] * ratio * np.absolute(ezs_sum)**2

#############################


rcs = open(directory+"/MRCS_exact_{}.txt".format(ratio), "w")
np.savetxt(rcs, MRCS)
rcs.close()


##rcs = open(directory+"/MRCS_exact_PEC.txt", "w")
##np.savetxt(rcs, MRCS)
##rcs.close()
################################################################################
