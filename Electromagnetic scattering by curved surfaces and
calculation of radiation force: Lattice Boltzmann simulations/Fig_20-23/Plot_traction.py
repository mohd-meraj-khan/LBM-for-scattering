import numpy as np ; from numpy.linalg import *
import sys
import os
import csv
import math
import scipy.special as sc
from math import e
import cmath
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec



a  = 100
velocity = 1.0/3

ratio91 = 0.91
wavelength91 = a / ratio91
period91 = wavelength91 / velocity

ratio93 = 0.93
wavelength93 = a / ratio93
period93 = wavelength93 / velocity

ratio95 = 0.95
wavelength95 = a / ratio95
period95 = wavelength95 / velocity



dphi = 0.1
phi = np.arange(0, 360 + dphi, dphi)

er = 4

tracExactAvg_091  = np.loadtxt('data/tractionExact_Avg_er_4_ratio_0.91.txt')
tracExactAvg_093  = np.loadtxt('data/tractionExact_Avg_er_4_ratio_0.93.txt')
tracExactAvg_095  = np.loadtxt('data/tractionExact_Avg_er_4_ratio_0.95.txt')

tracExactIns_091  = np.loadtxt('data/tractionExact_Ins_er_4_ratio_0.91.txt')
tracExactIns_093  = np.loadtxt('data/tractionExact_Ins_er_4_ratio_0.93.txt')
tracExactIns_095  = np.loadtxt('data/tractionExact_Ins_er_4_ratio_0.95.txt')




##tracLBM_091 = np.loadtxt('code/data/traction_ratio_0.91_a_100.txt')
##tracLBM_093 = np.loadtxt('code/data/traction_ratio_0.93_a_100.txt')
##tracLBM_095 = np.loadtxt('code/data/traction_ratio_0.95_a_100.txt')
##
##
##
##trac091 = open("data/tracLBM_091.txt", "w")
##np.savetxt(trac091, tracLBM_091[-int(period91):], fmt='%.4e')
##trac091.close()
##
##trac093 = open("data/tracLBM_093.txt", "w")
##np.savetxt(trac093, tracLBM_093[-int(period93):], fmt='%.4e')
##trac093.close()
##
##trac095 = open("data/tracLBM_095.txt", "w")
##np.savetxt(trac095, tracLBM_095[-int(period95):], fmt='%.4e')
##trac095.close()



tracLBM_091 = np.loadtxt('data/tracLBM_091.txt')
tracLBM_093 = np.loadtxt('data/tracLBM_093.txt')
tracLBM_095 = np.loadtxt('data/tracLBM_095.txt')

########################################################################################################################################################
########################################################################################################################################################
			
plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 1)
plt.rc('text', usetex = True)




##################################################################################################
##################################################################################################

fig = plt.figure(figsize = (6.69, 1.85), dpi=600)

gs=GridSpec(1,3)


ax0 = fig.add_subplot(gs[0,0])
ax0.plot(phi, tracExactAvg_091, 'k-')
ax0.plot(phi, np.average(tracLBM_091, axis=0), 'r--')


plt.text(30, 0.35, r'(a)', backgroundcolor='white')

plt.legend([r'Analytical', r'LBM'])

plt.xticks([0, 120, 240, 360])

plt.xlabel(r'$\phi (^o)$')
ax0.set_ylabel(r'$\frac{\langle f_x \rangle}{ \varepsilon_0 E_0^2}$')



#########

ax1 = fig.add_subplot(gs[0,1])
ax1.plot(phi, tracExactAvg_093, 'k-')
ax1.plot(phi, np.average(tracLBM_093, axis=0), 'r--')

plt.text(30, -2, r'(b)', backgroundcolor='white')

plt.xticks([0, 120, 240, 360])

plt.xlabel(r'$\phi (^o)$')
##ax.set_ylabel(r'$\frac{ F_x  / L}{\lambda \varepsilon_0 E_0^2}$')


########

ax2 = fig.add_subplot(gs[0,2])
ax2.plot(phi, tracExactAvg_095, 'k-')
ax2.plot(phi, np.average(tracLBM_095, axis=0), 'r--')

plt.text(30, -15, r'(c)', backgroundcolor='white')

plt.xticks([0, 120, 240, 360])

plt.xlabel(r'$\phi (^o)$')
##ax2.set_ylabel(r'$\frac{ F_x  / L}{\lambda \varepsilon_0 E_0^2}$')

plt.tight_layout()

plt.savefig('gallery/tracAvg_er_{}.svg'.format(er))
plt.close()
##################################################################################################
##################################################################################################



##################################################################################################
##################################################################################################

fig = plt.figure(figsize = (6.69, 1.85), dpi=600)

gs=GridSpec(1,3)


ax0 = fig.add_subplot(gs[0,0])
ax0.plot(phi, tracExactIns_091, 'k-')
ax0.plot(phi, tracLBM_091[27], 'r--')


plt.text(30, 0.5, r'(a)', backgroundcolor='white')

plt.legend([r'Analytical', r'LBM'])

plt.xticks([0, 120, 240, 360])

plt.xlabel(r'$\phi (^o)$')
ax0.set_ylabel(r'$\frac{ f_x}{ \varepsilon_0 E_0^2}$')



#########

ax1 = fig.add_subplot(gs[0,1])
ax1.plot(phi, tracExactIns_093, 'k-')
ax1.plot(phi, tracLBM_093[14], 'r--')

plt.text(30, -3.5, r'(b)', backgroundcolor='white')

plt.xticks([0, 120, 240, 360])

plt.xlabel(r'$\phi (^o)$')
##ax.set_ylabel(r'$\frac{ F_x  / L}{\lambda \varepsilon_0 E_0^2}$')


########

ax2 = fig.add_subplot(gs[0,2])
ax2.plot(phi, tracExactIns_095, 'k-')
ax2.plot(phi, tracLBM_095[11], 'r--')

plt.text(30, -20, r'(c)', backgroundcolor='white')

plt.xticks([0, 120, 240, 360])

plt.xlabel(r'$\phi (^o)$')
##ax2.set_ylabel(r'$\frac{ F_x  / L}{\lambda \varepsilon_0 E_0^2}$')

plt.tight_layout()

plt.savefig('gallery/tracIns_er_{}.svg'.format(er))
plt.close()
##################################################################################################
##################################################################################################

##sys.exit()

##################################################################################################
##################################################################################################

fig = plt.figure(figsize = (2.69, 1.85), dpi=600)
gs=GridSpec(1,1)
plt.ion()

for i in range(int(period91)):
    fig.clear()
    
    ax0 = fig.add_subplot(gs[0,0])
    ax0.plot(phi, tracExactIns_091, 'k-')
    ax0.plot(phi, tracLBM_091[i], 'r--')


##    plt.text(30, -1.75, r'(a)', backgroundcolor='white')

    plt.legend([r'Analytical', r'LBM'])

    plt.xticks([0, 120, 240, 360])

    plt.xlabel(r'$\phi (^o)$')
    ax0.set_ylabel(r'$\frac{ f_x  }{ \varepsilon_0 E_0^2}$')


    plt.tight_layout()

    plt.savefig('gallery/trac91_er_{}_{}.svg'.format(er, i))
plt.close()
##################################################################################################


##################################################################################################

fig = plt.figure(figsize = (2.69, 1.85), dpi=600)
gs=GridSpec(1,1)
plt.ion()

for i in range(int(period93)):
    fig.clear()
    
    ax0 = fig.add_subplot(gs[0,0])
    ax0.plot(phi, tracExactIns_093, 'k-')
    ax0.plot(phi, tracLBM_093[i], 'r--')


##    plt.text(30, -1.75, r'(a)', backgroundcolor='white')

    ##plt.legend([r'$ \langle $ Analytical $ \rangle $', r'LBM', r'$ \langle $ LBM $ \rangle $'])

    plt.xticks([0, 120, 240, 360])

    plt.xlabel(r'$\phi (^o)$')
##    ax0.set_ylabel(r'$\frac{ f_x  / L}{\lambda \varepsilon_0 E_0^2}$')


    plt.tight_layout()

    plt.savefig('gallery/trac93_er_{}_{}.svg'.format(er, i))
plt.close()
##################################################################################################



##################################################################################################

fig = plt.figure(figsize = (2.69, 1.85), dpi=600)
gs=GridSpec(1,1)
plt.ion()

for i in range(int(period95)):
    fig.clear()
    
    ax0 = fig.add_subplot(gs[0,0])
    ax0.plot(phi, tracExactIns_095, 'k-')
    ax0.plot(phi, tracLBM_095[i], 'r--')


##    plt.text(30, -1.75, r'(a)', backgroundcolor='white')

    ##plt.legend([r'$ \langle $ Analytical $ \rangle $', r'LBM', r'$ \langle $ LBM $ \rangle $'])

    plt.xticks([0, 120, 240, 360])

    plt.xlabel(r'$\phi (^o)$')
##    ax0.set_ylabel(r'$\frac{ f_x  / L}{\lambda \varepsilon_0 E_0^2}$')


    plt.tight_layout()

    plt.savefig('gallery/trac95_er_{}_{}.svg'.format(er, i))
plt.close()
##################################################################################################






##################################################################################################
##################################################################################################
##################################################################################################

fig = plt.figure(figsize = (2.69, 1.85), dpi=600)
gs=GridSpec(1,1)
plt.ion()


fig.clear()

ax0 = fig.add_subplot(gs[0,0])
ax0.plot(phi, tracExactAvg_091, 'k-')
ax0.plot(phi, np.average(tracLBM_091, axis=0), 'r--')


##plt.text(30, -1.75, r'(a)', backgroundcolor='white')

plt.legend([r'Analytical', r'LBM'])

plt.xticks([0, 120, 240, 360])

plt.xlabel(r'$\phi (^o)$')
ax0.set_ylabel(r'$\frac{ \langle f_x \rangle}{ \varepsilon_0 E_0^2}$')


plt.tight_layout()

plt.savefig('gallery/trac91_er_{}_avg.svg'.format(er))
plt.close()
##################################################################################################



##################################################################################################

fig = plt.figure(figsize = (2.69, 1.85), dpi=600)
gs=GridSpec(1,1)
plt.ion()


fig.clear()

ax0 = fig.add_subplot(gs[0,0])
ax0.plot(phi, tracExactAvg_091, 'k-')
ax0.plot(phi, np.average(tracLBM_091, axis=0), 'r--')


##plt.text(30, -1.75, r'(a)', backgroundcolor='white')

##plt.legend([r'$ \langle $ Analytical $ \rangle $', r'LBM', r'$ \langle $ LBM $ \rangle $'])

plt.xticks([0, 120, 240, 360])

plt.xlabel(r'$\phi (^o)$')
##ax0.set_ylabel(r'$\frac{ f_x  / L}{\lambda \varepsilon_0 E_0^2}$')


plt.tight_layout()

plt.savefig('gallery/trac93_er_{}_avg.svg'.format(er))
plt.close()
##################################################################################################



##################################################################################################

fig = plt.figure(figsize = (2.69, 1.85), dpi=600)
gs=GridSpec(1,1)
plt.ion()


fig.clear()

ax0 = fig.add_subplot(gs[0,0])
ax0.plot(phi, tracExactAvg_095, 'k-')
ax0.plot(phi, np.average(tracLBM_095, axis=0), 'r--')


##plt.text(30, -1.75, r'(a)', backgroundcolor='white')

##plt.legend([r'$ \langle $ Analytical $ \rangle $', r'LBM', r'$ \langle $ LBM $ \rangle $'])

plt.xticks([0, 120, 240, 360])

plt.xlabel(r'$\phi (^o)$')
##ax0.set_ylabel(r'$\frac{ f_x  / L}{\lambda \varepsilon_0 E_0^2}$')


plt.tight_layout()

plt.savefig('gallery/trac95_er_{}_avg.svg'.format(er))
plt.close()
##################################################################################################








