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



er2 = 4

r = 5

a = 1


dratio = 0.001
ratio = np.arange(0.9, 1.1 + dratio, dratio)

dphi = 1
phi = np.arange(0, 360 + 0.001, dphi)

BRCS = np.zeros(len(phi))


er1 = 1
mur1 = 1
mur2 = 1

E0 = 1

v1 = 1 / (3 * np.sqrt(er1 * mur1))
v2 = 1 / (3 * np.sqrt(er2 * mur2))

alpha = np.sqrt((er1 * mur2) / (er2 * mur1))


C_sca = np.zeros(len(ratio))
Q_sca = np.zeros(len(ratio))

########################################## TSCS ###################################

for m in range(len(ratio)):

    k1a = 2 * np.pi * ratio[m]
    k2a = 2 * np.pi * v1 / v2 * ratio[m]
    k1r = 2 * np.pi * ratio[m] * r

    w_sca = np.zeros(len(phi))
    
    n = 1000

    for i in range(len(phi)):
       
        ez_s_sum     = 0
        hr_s_sum     = 0
        hphi_s_sum   = 0
            
        for l in range(0, n+1, 1):
                
            num  = (alpha * sc.jvp(l, k1a) * sc.jv(l, k2a) - sc.jv(l, k1a) * sc.jvp(l, k2a))
            den  = (sc.hankel2(l, k1a) * sc.jvp(l, k2a) - alpha * sc.h2vp(l, k1a) * sc.jv(l, k2a))
                
##            num = - sc.jv(l, k1a)
##            den = sc.hankel2(l, k1a)

            Cn = num / den

            if (l == 0):
                ez_s   = (0 - 1j)**l * (Cn * sc.hankel2(l, k1r)) * np.cos(l*phi[i]*np.pi/180)
                hr_s   = (0 - 1j)**l * (Cn * sc.hankel2(l, k1r)) * np.cos(l*phi[i]*np.pi/180) * l
                hphi_s = (0 - 1j)**l * (Cn * sc.h2vp(l, k1r)) * np.cos(l*phi[i]*np.pi/180)
            else :
                ez_s   = 2 * (0 - 1j)**l * (Cn * sc.hankel2(l, k1r)) * np.cos(l*phi[i]*np.pi/180)
                hr_s   = 2 * (0 - 1j)**l * (Cn * sc.hankel2(l, k1r)) * np.cos(l*phi[i]*np.pi/180) * l
                hphi_s = 2 * (0 - 1j)**l * (Cn * sc.h2vp(l, k1r)) * np.cos(l*phi[i]*np.pi/180)
                
            ez_s_sum     += ez_s
            hr_s_sum     += hr_s
            hphi_s_sum   += hphi_s 

            if (np.absolute(hphi_s) <= 1e-10 and np.absolute(np.cos(l*phi[i]*np.pi/180)) >= 1e-10):
##                print(l)
                break       

        Ez_s   = ez_s_sum
        Hr_s   = hr_s_sum * (-1)  / (2 * np.pi * mur1 * ratio[m] * r)
        Hphi_s = hphi_s_sum * (0 - 1j) * np.sqrt(er1 / mur1)
        
        w_sca[i] = - 0.5 * np.real(Ez_s * np.conjugate(Hphi_s))


    C_sca[m] = np.sum(w_sca) * r*a * (dphi*np.pi/180) / (0.5 * np.sqrt(er1 / mur1) * E0**2)
    Q_sca[m] = C_sca[m] / (2 * a)


print(Q_sca)

###################################

total_scs = open(directory+"/Q_sca_exact_{}.txt".format(er2), "w")
np.savetxt(total_scs, Q_sca)
total_scs.close()

###################################################################################



