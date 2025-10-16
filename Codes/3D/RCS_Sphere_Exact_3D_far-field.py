import numpy as np
import matplotlib.pyplot as plt
import sys
import os

from scipy.interpolate import RectBivariateSpline

import math
import scipy.special as sc
from math import e
import cmath
from scipy.optimize import curve_fit

from Module_Parameters_3D import *

directory_rcs = 'data/rcs'
if not os.path.exists(directory_rcs):
    os.makedirs(directory_rcs)



'''see equation 10.31 of Balanis for spherical Bessel and Hankel functions'''



'''functions for spherical Bessel and Hankel functions'''
def Jn_x(n, x):
    result = np.sqrt(np.pi*x / (2)) * sc.jv(n+0.5, x)
    return result

def Hn2_x(n, x):
    result = np.sqrt(np.pi*x / (2)) * sc.hankel2(n+0.5, x)
    return result


'''functions for first derivatives of spherical Bessel and Hankel functions'''
def dJn_dx(n, x):
    result = np.sqrt(np.pi/(8*x)) * (x * sc.jv(n-0.5, x) + sc.jv(n+0.5, x) - x * sc.jv(n+1.5, x))
    return result

def dHn2_dx(n, x):
    result = np.sqrt(np.pi/(8*x)) * (x * sc.hankel2(n-0.5, x) + sc.hankel2(n+0.5, x) - x * sc.hankel2(n+1.5, x))
    return result


'''functions for second derivatives of spherical Bessel and Hankel functions'''
def d2Jn_dx2(n, x):
    result = np.sqrt(np.pi/(32*x**3)) * (x**2 * sc.jv(n-1.5, x) + 2*x * sc.jv(n-0.5, x) - sc.jv(n+0.5, x) - 2*x**2 * sc.jv(n+0.5, x) - 2*x * sc.jv(n+1.5, x) + x**2 * sc.jv(n+2.5, x))
    return result

def d2Hn2_dx2(n, x):
    result = np.sqrt(np.pi/(32*x**3)) * (x**2 * sc.hankel2(n-1.5, x) + 2*x * sc.hankel2(n-0.5, x) - sc.hankel2(n+0.5, x) - 2*x**2 * sc.hankel2(n+0.5, x) - 2*x * sc.hankel2(n+1.5, x) + x**2 * sc.hankel2(n+2.5, x))
    return result


'''functions of associated Legendre function and derivative'''
def Pn_1(n, x):
    P = np.array(sc.lpmn(1, n, x))
    return P[0, 1, n]

def dPn_1(n, x):
    P = np.array(sc.lpmn(1, n, x))
    return P[1, 1, n]




##########################################################################################################################################
######                                            EXACT SOLUTION OF SCATTERED ELECTRIC FIRLD                                          ####
##########################################################################################################################################

epr = er2
Ratio = ratio

mur = 1

N = 1000

a = 50

d2r = np.pi / 180
Theta = np.arange(0, 360+1, 1)



A_theta = np.zeros(len(Theta), dtype=complex)
A_phi   = np.zeros(len(Theta), dtype=complex)

BRCS      = np.zeros(len(Theta))


E0 = 1



wavelength = a / Ratio
k0 = 2 * np.pi / wavelength

kd = k0 * np.sqrt(epr*mur)
     



phi = phi0

for i in range(len(Theta)):
    theta = Theta[i]

    sum_a_theta = 0
    sum_a_phi   = 0
                
    for n in range(1,N+1,1):

        a_theta = 0
        a_phi   = 0


        if (Theta[i] == 180):
            An = (-1)**n * n*(n+1)/2
            Bn = (-1)**n * n*(n+1)/2
        elif (Theta[i] == 0 or Theta[i] == 360):
            An = - n*(n+1)/2
            Bn =   n*(n+1)/2
        else:
            An = Pn_1(n, np.cos(theta*d2r)) / np.sin(theta*d2r)
            Bn = np.sin(theta*d2r) * dPn_1(n, np.cos(theta*d2r))

                    
        an = (0 + 1j)**(-n) * (2*n + 1) / (n * (n + 1))

        # coefficients for PEC sphere
        bn = -an * dJn_dx(n, k0*a) / dHn2_dx(n, k0*a)
        cn = -an * Jn_x(n, k0*a) / Hn2_x(n, k0*a)

        # coefficients for dielectric sphere
        bn = an * (- np.sqrt(epr) * dJn_dx(n, k0*a) * Jn_x(n, kd*a) + np.sqrt(mur) * Jn_x(n, k0*a) * dJn_dx(n, kd*a)
                    ) / (np.sqrt(epr) * dHn2_dx(n, k0*a) * Jn_x(n, kd*a) - np.sqrt(mur) * Hn2_x(n, k0*a) * dJn_dx(n, kd*a))
                        
        cn = an * (- np.sqrt(epr) * Jn_x(n, k0*a) * dJn_dx(n, kd*a) + np.sqrt(mur) * dJn_dx(n, k0*a) * Jn_x(n, kd*a)
                    ) / (np.sqrt(epr) * Hn2_x(n, k0*a) * dJn_dx(n, kd*a) - np.sqrt(mur) * dHn2_dx(n, k0*a) * Jn_x(n, kd*a))


        a_theta = (0 + 1j)**n * (bn * Bn - cn * An )
        a_phi   = (0 + 1j)**n * (bn * An - cn * Bn )

        sum_a_theta += a_theta
        sum_a_phi   += a_phi
                    
        cutoff = 1e-10

        if (np.absolute(a_theta) <= cutoff and np.absolute(a_phi) <= cutoff):
            break

    
    A_theta[i] = sum_a_theta
    A_phi[i]   = sum_a_phi


    BRCS[i] = 1 / (np.pi**2 * Ratio**2) * ( np.cos(phi*d2r)**2 * np.absolute(A_theta[i])**2 + np.sin(phi*d2r)**2 * np.absolute(A_phi[i])**2 )


np.save(directory_rcs+'/BRCS_exact_{}_{}_{}.npy'.format(Ratio, epr, phi), BRCS)

##########################################################################################################################################




##################################
plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)
##################################


############################################# Radial snapshot #############################################

##fig, ax = plt.subplots(figsize = (3.35, 2.05), dpi=600, constrained_layout = True)
##
##plt.plot(Theta, BRCS,'k-')
##
####plt.plot(Theta, 10*np.log10(BRCS),'k')
##
##plt.xlabel(r'$\theta (^o)$')
##plt.ylabel(r'$ \sigma /(\pi a^2)$')
##
##plt.yscale('log')
##
##    
##plt.grid(True)
##ax.grid(which='minor')
##    
##ax.set_xticks([0, 90, 180, 270, 360])
##ax.tick_params(which = 'both', direction='in')
##
##plt.savefig('BRCS_Exact_{}_{}_{}.svg'.format(Ratio, epr, phi))
##plt.close(fig)
    ##########################################################################################################


sys.exit()



##########################################################################################################################################
######                                            EXACT SOLUTION OF SCATTERED ELECTRIC FIRLD                                          ####
##########################################################################################################################################

epr = 2
mur = 1

N = 1000

a = 50

d2r = np.pi / 180

Ratio = np.arange(0.02, 1.6, 0.01)

Theta = 180
phi   = 180



A_theta = np.zeros(len(Ratio), dtype=complex)
A_phi   = np.zeros(len(Ratio), dtype=complex)

MRCS      = np.zeros(len(Ratio))


E0 = 1


for i in range(len(Ratio)):
    ratio = Ratio[i]

    wavelength = a / ratio
    k0 = 2 * np.pi / wavelength
    kd = k0 * np.sqrt(epr*mur)

    
    sum_a_theta = 0
    sum_a_phi   = 0
                
    for n in range(1,N+1,1):

        a_theta = 0
        a_phi   = 0


        if (Theta == 180):
            An = (-1)**n * n*(n+1)/2
            Bn = (-1)**n * n*(n+1)/2
        elif (Theta == 0 or Theta == 360):
            An = - n*(n+1)/2
            Bn =   n*(n+1)/2
        else:
            An = Pn_1(n, np.cos(theta*d2r)) / np.sin(theta*d2r)
            Bn = np.sin(theta*d2r) * dPn_1(n, np.cos(theta*d2r))

                    
        an = (0 + 1j)**(-n) * (2*n + 1) / (n * (n + 1))

        # coefficients for PEC sphere
        bn = -an * dJn_dx(n, k0*a) / dHn2_dx(n, k0*a)
        cn = -an * Jn_x(n, k0*a) / Hn2_x(n, k0*a)

        # coefficients for dielectric sphere
##        bn = an * (- np.sqrt(epr) * dJn_dx(n, k0*a) * Jn_x(n, kd*a) + np.sqrt(mur) * Jn_x(n, k0*a) * dJn_dx(n, kd*a)
##                    ) / (np.sqrt(epr) * dHn2_dx(n, k0*a) * Jn_x(n, kd*a) - np.sqrt(mur) * Hn2_x(n, k0*a) * dJn_dx(n, kd*a))
##                        
##        cn = an * (- np.sqrt(epr) * Jn_x(n, k0*a) * dJn_dx(n, kd*a) + np.sqrt(mur) * dJn_dx(n, k0*a) * Jn_x(n, kd*a)
##                    ) / (np.sqrt(epr) * Hn2_x(n, k0*a) * dJn_dx(n, kd*a) - np.sqrt(mur) * dHn2_dx(n, k0*a) * Jn_x(n, kd*a))


        a_theta = (0 + 1j)**n * (bn * Bn - cn * An )
        a_phi   = (0 + 1j)**n * (bn * An - cn * Bn )



        sum_a_theta += a_theta
        sum_a_phi   += a_phi
                    
        cutoff = 1e-10
                    
        if (np.absolute(a_theta) <= cutoff and np.absolute(a_phi) ):
##            print(n)
            break

    A_theta[i] = sum_a_theta
    A_phi[i]   = sum_a_phi


    MRCS[i] = 1 / (np.pi**2 * ratio**2) * ( np.cos(phi*d2r)**2 * np.absolute(A_theta[i])**2 + np.sin(phi*d2r)**2 * np.absolute(A_phi[i])**2 )



np.save(directory_rcs+'/MRCS_exact_{}.npy'.format(epr), MRCS)

##########################################################################################################################################





plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)


############################################# Radial snapshot #############################################

fig, ax = plt.subplots(figsize = (3.35, 1.75), dpi=600, constrained_layout = True)

plt.plot(Ratio, MRCS,'k')

##plt.plot(Theta, 10*np.log10(MRCS),'k')


plt.xlabel(r'$a / \lambda$')
plt.ylabel(r'$ \sigma /(\pi a^2)$')

plt.yscale('log')
##plt.xscale('log')

    
plt.grid(True)
ax.grid(which='minor')
    
ax.tick_params(which = 'both', direction='in')

plt.savefig('MRCS_Exact_{}.svg'.format(epr))
plt.close(fig)
##########################################################################################################


