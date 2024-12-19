import numpy as np ; from numpy.linalg import *
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.patches as patches
import sys
import os
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.gridspec import GridSpec
from matplotlib.patches import ConnectionPatch
import csv
import matplotlib as mpl




##N = [1, 2, 3, 4, 5]
##epsilon = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
##
##

a = 200

Nx = 8*a
Ny = 8*a





########### Plot parameters ################
plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)
###########################################
    
fig, axs = plt.subplots(1, 3, figsize = (6.69, 2.45), constrained_layout=True, dpi=600)

plt.ion()




fig.suptitle(r'$H_z^{tot}$')


nx = Nx//(4*a)
ny = Ny//(4*a)


A05 = np.loadtxt("Hz_0.5.txt")
A10 = np.loadtxt("Hz_1.txt")
A15 = np.loadtxt("Hz_1.5.txt")


im1 = axs[0].imshow(A05, vmin = -2, vmax = 2, cmap='seismic', origin='lower')

axs[0].set_xticks(np.linspace(0,Nx,5))
axs[0].set_xticklabels([0, nx, 2*nx, 3*nx, 4*nx])
axs[0].set_xlabel(r'$x/a$')
        
axs[0].set_yticks(np.linspace(0,Ny,5))
axs[0].set_yticklabels([0, ny, 2*ny, 3*ny, 4*ny])
axs[0].set_ylabel(r'$y/a$')

axs[0].text(0.05, 0.9, r'$A = 0.5$', fontsize=10, transform=axs[0].transAxes, bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=0.1'))


im2 = axs[1].imshow(A10, vmin = -2, vmax = 2, cmap='seismic', origin='lower')

axs[1].set_xticks(np.linspace(0,Nx,5))
axs[1].set_xticklabels([0, nx, 2*nx, 3*nx, 4*nx])
axs[1].set_xlabel(r'$x/a$')
        
axs[1].set_yticks(np.linspace(0,Ny,5))
axs[1].set_yticklabels([])
##axs[1].set_ylabel(r'$x$')

axs[1].text(0.05, 0.9, r'$A = 1.0$', fontsize=10, transform=axs[1].transAxes, bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=0.1'))


im3 = axs[2].imshow(A15, vmin = -2, vmax = 2, cmap='seismic', origin='lower')

axs[2].set_xticks(np.linspace(0,Nx,5))
axs[2].set_xticklabels([0, nx, 2*nx, 3*nx, 4*nx])
axs[2].set_xlabel(r'$x/a$')
        
axs[2].set_yticks(np.linspace(0,Ny,5))
axs[2].set_yticklabels([])
##axs[2].set_ylabel(r'$x$')

axs[2].text(0.05, 0.9, r'$A = 1.5$', fontsize=10, transform=axs[2].transAxes, bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=0.1'))


fig.colorbar(im3, ax = axs[2], location='right', shrink=1.0, aspect=25)


##axs[j].text(0.05, 0.9, r'$N = {}$'.format(N[j]), fontsize=8, transform=axs[j].transAxes, bbox=dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.3'))
        
plt.savefig("corrugatedSnap.svg")
plt.close()




########### Plot parameters ################
plt.rc('font', family = 'serif', size = 20)
plt.rc('xtick', labelsize = 20)
plt.rc('ytick', labelsize = 20)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)
###########################################


################################################################

fig, axs = plt.subplots(1, 1, figsize = (8.0139, 6.2739), constrained_layout=True, dpi=600)

plt.ion()




fig.suptitle(r'$H_z^{tot}$')


nx = Nx//(4*a)
ny = Ny//(4*a)


A05 = np.loadtxt("Hz_0.5.txt")



im1 = axs.imshow(A05, vmin = -2, vmax = 2, cmap='seismic', origin='lower')

axs.set_xticks(np.linspace(0,Nx,5))
axs.set_xticklabels([0, nx, 2*nx, 3*nx, 4*nx])
axs.set_xlabel(r'$x/a$')
        
axs.set_yticks(np.linspace(0,Ny,5))
axs.set_yticklabels([0, ny, 2*ny, 3*ny, 4*ny])
axs.set_ylabel(r'$y/a$')

axs.text(0.05, 0.9, r'$A = 0.5$', transform=axs.transAxes, bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=0.1'))



fig.colorbar(im1, ax = axs, location='right', shrink=1.0, aspect=25)


plt.savefig("Highlight_Corrugated_0.5.eps")
plt.close()
                

################################################################

fig, axs = plt.subplots(1, 1, figsize = (8.0139, 6.2739), constrained_layout=True, dpi=600)

plt.ion()




fig.suptitle(r'$H_z^{tot}$')


nx = Nx//(4*a)
ny = Ny//(4*a)


A10 = np.loadtxt("Hz_1.txt")


im1 = axs.imshow(A10, vmin = -2, vmax = 2, cmap='seismic', origin='lower')

axs.set_xticks(np.linspace(0,Nx,5))
axs.set_xticklabels([0, nx, 2*nx, 3*nx, 4*nx])
axs.set_xlabel(r'$x/a$')
        
axs.set_yticks(np.linspace(0,Ny,5))
axs.set_yticklabels([0, ny, 2*ny, 3*ny, 4*ny])
axs.set_ylabel(r'$y/a$')

axs.text(0.05, 0.9, r'$A = 1.0$', transform=axs.transAxes, bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=0.1'))


fig.colorbar(im1, ax = axs, location='right', shrink=1.0, aspect=25)

        
plt.savefig("Highlight_Corrugated_1.0.eps")
plt.close()


################################################################

fig, axs = plt.subplots(1, 1, figsize = (8.0139, 6.2739), constrained_layout=True, dpi=600)

plt.ion()




fig.suptitle(r'$H_z^{tot}$')


nx = Nx//(4*a)
ny = Ny//(4*a)


A15 = np.loadtxt("Hz_1.5.txt")


im1 = axs.imshow(A15, vmin = -2, vmax = 2, cmap='seismic', origin='lower')

axs.set_xticks(np.linspace(0,Nx,5))
axs.set_xticklabels([0, nx, 2*nx, 3*nx, 4*nx])
axs.set_xlabel(r'$x/a$')
        
axs.set_yticks(np.linspace(0,Ny,5))
axs.set_yticklabels([0, ny, 2*ny, 3*ny, 4*ny])
axs.set_ylabel(r'$y/a$')

axs.text(0.05, 0.9, r'$A = 1.5$', transform=axs.transAxes, bbox=dict(facecolor='white', edgecolor='none', boxstyle='square,pad=0.1'))



fig.colorbar(im1, ax = axs, location='right', shrink=1.0, aspect=25)


plt.savefig("Highlight_Corrugated_1.5.eps")
plt.close()
                
