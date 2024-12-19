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

epr = 2
mur = 1

N = 1000

a = 50

d2r = np.pi / 180

Ratio = [0.5]#np.arange(0.03, 2, 0.001)
R     = [2]
Theta = np.arange(0, 180+1, 1)
Phi   = [0] 


Es_r     = np.zeros((len(Ratio), len(R), len(Theta), len(Phi)), dtype=complex)
Es_theta = np.zeros((len(Ratio), len(R), len(Theta), len(Phi)), dtype=complex)
Es_phi   = np.zeros((len(Ratio), len(R), len(Theta), len(Phi)), dtype=complex)

RCS      = np.zeros((len(Ratio), len(R), len(Theta), len(Phi)))


E0 = 1

for s in range(len(Ratio)):

    wavelength = a / Ratio[s]
    k0 = 2 * np.pi / wavelength

    kd = k0 * np.sqrt(epr*mur)
    
    for i in range(len(R)):
        r = R[i]*a
        for m in range(len(Theta)):
            theta = Theta[m]
            for l in range(len(Phi)):
                phi = Phi[l]
            
                sum_es_r     = 0
                sum_es_theta = 0
                sum_es_phi   = 0
                
                for n in range(1,N+1,1):

                    es_r     = 0
                    es_theta = 0
                    es_phi   = 0


                    if (Theta[m] == 180):
                        An = (-1)**n * n*(n+1)/2
                        Bn = (-1)**n * n*(n+1)/2
                    elif (Theta[m] == 0):
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
##                    bn = an * (- np.sqrt(epr) * dJn_dx(n, k0*a) * Jn_x(n, kd*a) + np.sqrt(mur) * Jn_x(n, k0*a) * dJn_dx(n, kd*a)
##                               ) / (np.sqrt(epr) * dHn2_dx(n, k0*a) * Jn_x(n, kd*a) - np.sqrt(mur) * Hn2_x(n, k0*a) * dJn_dx(n, kd*a))
##                        
##                    cn = an * (- np.sqrt(epr) * Jn_x(n, k0*a) * dJn_dx(n, kd*a) + np.sqrt(mur) * dJn_dx(n, k0*a) * Jn_x(n, kd*a)
##                               ) / (np.sqrt(epr) * Hn2_x(n, k0*a) * dJn_dx(n, kd*a) - np.sqrt(mur) * dHn2_dx(n, k0*a) * Jn_x(n, kd*a))




                    es_r     = bn * (d2Hn2_dx2(n, k0*r) + Hn2_x(n, k0*r)) * Pn_1(n, np.cos(theta*d2r))
                    es_theta = (0 + 1j) * bn * dHn2_dx(n, k0*r) * Bn - cn * Hn2_x(n, k0*r) * An
                    es_phi   = (0 + 1j) * bn * dHn2_dx(n, k0*r) * An - cn * Hn2_x(n, k0*r) * Bn

                    sum_es_r     += es_r
                    sum_es_theta += es_theta
                    sum_es_phi   += es_phi
                    
                    cutoff = 1e-10
                    
                    if (np.absolute(es_r) <= cutoff and np.absolute(es_theta) <= cutoff and np.absolute(es_phi) <= cutoff):
##                        print(n)
                        break

                Es_r[s, i, m, l]     = (0 - 1j) * E0 * np.cos(phi*d2r) * sum_es_r
                Es_theta[s, i, m, l] = E0 / (k0*r) * np.cos(phi*d2r) * sum_es_theta
                Es_phi[s, i, m, l]   = E0 / (k0*r) * np.sin(phi*d2r) * sum_es_phi



                RCS[s, i, m, l] = 4*r**2/a**2 * (np.absolute(Es_r[s, i, m, l])**2 + np.absolute(Es_theta[s, i, m, l])**2 + np.absolute(Es_phi[s, i, m, l])**2)

Es_r     = np.absolute(Es_r)
Es_theta = np.absolute(Es_theta)
Es_phi   = np.absolute(Es_phi)



brcs = open(directory+"/BRCS_Exact_er_{}.txt".format(epr), "w")
np.savetxt(brcs, RCS[0, 0, :, 0])
brcs.close()


brcs = open(directory+"/BRCS_Exact_er_{}_{}.txt".format(epr, R[0]), "w")
np.savetxt(brcs, RCS[0, 0, :, 0])
brcs.close()

##########################################################################################################################################





plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)


############################################# Radial snapshot #############################################

fig, ax = plt.subplots(figsize = (3.35, 2.25), dpi=600, constrained_layout = True)

plt.plot(Theta, RCS[0, 0, :, 0],'k')

##plt.plot(Theta, 10*np.log10(RCS[0, 0, :, 0]),'k')


plt.xlabel(r'$a / \lambda$')
plt.ylabel(r'$ \sigma /(\pi a^2)$')

plt.yscale('log')
    
plt.legend(['my plot','Balanis'])

    
plt.grid(True)
ax.grid(which='minor')
    
##plt.title(r'PEC sphere, $r/a = {}$, $k_0 a = {}$, $\phi = {}$'.format(R[i], k0a, phi))
##plt.xticks(np.linspace(0,180,7))
ax.tick_params(which = 'both', direction='in')
##fig.tight_layout()
plt.savefig('RCS_sphere_er_{}.svg'.format(epr))
plt.close(fig)
##########################################################################################################



