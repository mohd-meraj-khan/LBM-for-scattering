import numpy as np
import matplotlib.pyplot as plt


import math
import scipy.special as sc
from math import e
import cmath





ratio = 1

a=100


dphi = 0.1
phi = np.arange(0, 360 + 0.01, dphi)




#################### Material properties ##############
er1 = 1        # value of the permittivity in domain 1
mur1 = 1       # value of the permeability in domain 1

er2 = 4        # value of the permittivity in domain 2
mur2 = 1       # value of the permeability in domain 2
#######################################################

################ Wavelength parameters ################
v1 = 1 / (3 * np.sqrt(er1 * mur1))
v2 = 1 / (3 * np.sqrt(er2 * mur2))

omega = 2 * np.pi * ratio  * v1 / a
period = int(np.round(a / (v1 * ratio)))

alpha = np.sqrt((er1 * mur2) / (er2 * mur1))

n = 1000


k1a = 2 * np.pi * ratio
k2a = 2 * np.pi * v1 / v2 * ratio


print(period)

############################################################################################################################################################################
########################################################## Instantanuous traction vector in x directionn ###################################################################

##t = [0, int(np.round(period/8)), int(np.round(period/4)), int(np.round(period/2))]
t = [50]#[int(np.round(period/4))]

EzTot   = np.zeros((len(phi), len(t)))
HrTot   = np.zeros((len(phi), len(t)))
HphiTot = np.zeros((len(phi), len(t)))

fr   = np.zeros((len(phi), len(t)))
fphi = np.zeros((len(phi), len(t)))
fx   = np.zeros((len(phi), len(t)))


for i in range(len(phi)):
    eztot_sum   = 0
    hrtot_sum   = 0
    hphitot_sum = 0
    
    for l in range(0, n+1, 1):
        
        num  = (alpha * sc.jvp(l, k1a) * sc.jv(l, k2a) - sc.jv(l, k1a) * sc.jvp(l, k2a))
        den  = (sc.hankel2(l, k1a) * sc.jvp(l, k2a) - alpha * sc.h2vp(l, k1a) * sc.jv(l, k2a))

##                num = - sc.jv(l, k2a)
##                den = sc.hankel2(l, k2a)

        Cn = num / den

        if (l == 0):
            eztot   = (0 - 1j)**l * (sc.jv(l, k1a) + Cn * sc.hankel2(l, k1a)) * np.cos(l*phi[i]*np.pi/180)
            hrtot   = (0 - 1j) / (3*mur1*omega*a) * l * (0 - 1j)**l * (sc.jv(l, k1a) + Cn * sc.hankel2(l, k1a)) * np.sin(l*phi[i]*np.pi/180)
            hphitot = (0 - 1j) / (3*mur1*v1) * (0 - 1j)**l * (sc.jvp(l, k1a) + Cn * sc.h2vp(l, k1a)) * np.cos(l*phi[i]*np.pi/180)
        else :
            eztot   = 2 * (0 - 1j)**l * (sc.jv(l, k1a) + Cn * sc.hankel2(l, k1a)) * np.cos(l*phi[i]*np.pi/180)
            hrtot   = 2 * (0 - 1j) / (3*mur1*omega*a) * l * (0 - 1j)**l * (sc.jv(l, k1a) + Cn * sc.hankel2(l, k1a)) * np.sin(l*phi[i]*np.pi/180)
            hphitot = 2 * (0 - 1j) / (3*mur1*v1) * (0 - 1j)**l * (sc.jvp(l, k1a) + Cn * sc.h2vp(l, k1a)) * np.cos(l*phi[i]*np.pi/180)

        eztot_sum   += eztot
        hrtot_sum   += hrtot
        hphitot_sum += hphitot

        if (np.absolute(hphitot) <= 1e-10 and np.absolute(np.cos(l*phi[i]*np.pi/180)) >= 1e-10):
##                    print(l)
            break
        
    for g in range(len(t)):
        EzTot[i, g]   = np.real(eztot_sum * np.exp((0 + 1j)*omega * t[g]))
        HrTot[i, g]   = np.real(hrtot_sum * np.exp((0 + 1j)*omega * t[g]))
        HphiTot[i, g] = np.real(hphitot_sum * np.exp((0 + 1j)*omega * t[g]))
        

        fr[i, g]   = - 0.5 * EzTot[i, g]**2 + 0.5 * (HrTot[i, g]**2 - HphiTot[i, g]**2)
        fphi[i, g] = HphiTot[i, g] * HrTot[i, g]


        fx[i, g] = (fr[i, g] * np.cos(phi[i] * np.pi / 180) - fphi[i, g] * np.sin(phi[i] * np.pi / 180))

        

###################################


force = open("tractionExact_Ins_er_{}_ratio_{}.txt".format(er2, ratio), "w")
np.savetxt(force, fx[:, 0])
force.close()
############################################################################################################################################################################













############################################################################################################################################################################
################################################# Traction vector in x directionn averaged over one time period ############################################################

frAvg   = np.zeros(len(phi))
fphiAvg = np.zeros(len(phi))
fxAvg   = np.zeros(len(phi))

for i in range(len(phi)):
    eztot_sum   = 0
    hrtot_sum   = 0
    hphitot_sum = 0
    
    for l in range(0, n+1, 1):
        
        num  = (alpha * sc.jvp(l, k1a) * sc.jv(l, k2a) - sc.jv(l, k1a) * sc.jvp(l, k2a))
        den  = (sc.hankel2(l, k1a) * sc.jvp(l, k2a) - alpha * sc.h2vp(l, k1a) * sc.jv(l, k2a))

##        num = - sc.jv(l, k2a)
##        den = sc.hankel2(l, k2a)

        Cn = num / den

        if (l == 0):
            eztot   = (0 - 1j)**l * (sc.jv(l, k1a) + Cn * sc.hankel2(l, k1a)) * np.cos(l*phi[i]*np.pi/180)
            hrtot   = (0 - 1j) / (3*mur1*omega*a) * l * (0 - 1j)**l * (sc.jv(l, k1a) + Cn * sc.hankel2(l, k1a)) * np.sin(l*phi[i]*np.pi/180)
            hphitot = (0 - 1j) / (3*mur1*v1) * (0 - 1j)**l * (sc.jvp(l, k1a) + Cn * sc.h2vp(l, k1a)) * np.cos(l*phi[i]*np.pi/180)
        else :
            eztot   = 2 * (0 - 1j)**l * (sc.jv(l, k1a) + Cn * sc.hankel2(l, k1a)) * np.cos(l*phi[i]*np.pi/180)
            hrtot   = 2 * (0 - 1j) / (3*mur1*omega*a) * l * (0 - 1j)**l * (sc.jv(l, k1a) + Cn * sc.hankel2(l, k1a)) * np.sin(l*phi[i]*np.pi/180)
            hphitot = 2 * (0 - 1j) / (3*mur1*v1) * (0 - 1j)**l * (sc.jvp(l, k1a) + Cn * sc.h2vp(l, k1a)) * np.cos(l*phi[i]*np.pi/180)

        eztot_sum   += eztot
        hrtot_sum   += hrtot
        hphitot_sum += hphitot

        if (np.absolute(hphitot) <= 1e-10 and np.absolute(np.cos(l*phi[i]*np.pi/180)) >= 1e-10):
##            print(l)
            break
    

    frAvg[i]   = 0.5 * np.real( - 0.5 * np.absolute(eztot_sum)**2 + 0.5 * (np.absolute(hrtot_sum)**2 - np.absolute(hphitot_sum)**2))
    fphiAvg[i] = 0.5 * np.real(hphitot_sum * np.conjugate(hrtot_sum))


    fxAvg[i] = (frAvg[i] * np.cos(phi[i] * np.pi / 180) - fphiAvg[i] * np.sin(phi[i] * np.pi / 180))

    
    
###################################

FxSum = 0

for i in range(len(phi)):
    FxSum += fxAvg[i] * a * (dphi * np.pi / 180) / (a / ratio)

print(FxSum)

###################################


plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)


fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2), dpi=600, constrained_layout = True)

plt.plot(phi, fx[:,0], 'k-', phi, fxAvg, 'k--')

plt.xticks([0, 60, 120, 180,240, 300, 360])

plt.xlabel(r'$\phi (^o)$')
plt.ylabel(r'$\frac{f_x}{\varepsilon_0 E_0^2}, \frac{\left< f_x \right>}{\varepsilon_0 E_0^2}$')

plt.legend([r'Instantaneous', r'Averaged over 1 T'])
    
plt.savefig('tractionExact_er_{}_ratio_{}.svg'.format(er2, ratio))
plt.close()
############################################################################################################################################################################



force = open("tractionExact_Avg_er_{}_ratio_{}.txt".format(er2, ratio), "w")
np.savetxt(force, fxAvg)
force.close()
