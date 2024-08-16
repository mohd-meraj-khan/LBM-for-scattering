from ctypes import CDLL, POINTER
from ctypes import c_size_t, c_double, c_float, c_int
import numpy as np
from scipy.interpolate import RectBivariateSpline
import sys
import os
import time

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Wedge
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import ConnectionPatch

from Module_Traction import tractionFx, tractionFy, tractionFz
from Module_EM_Wave import planeWaveTM
import Module_Geometry

from parameters import *
from Shared_Lib import *


t0 = time.time()


pictures = 'gallery'
if not os.path.exists(pictures):
    os.makedirs(pictures)




######################
# LBM properties (DO NOT CHANGE)
Q = 7               # number of velocities at a grid
velocity = 1.0/3    # velocity of EM wave in vacuum
######################


######### DON'T CHANGE ########
wavelength = a / ratio
period = wavelength / velocity
omega = 2 * np.pi / period
###############################



Time = 3 * (Nx * np.sqrt(er1) + noOfReflections * 2 * a * np.sqrt(er2)) + noOfPeriods * period

print("Number of time steps :", int(Time))




###############################################################################################################
####      DEFINING AND INITILIZING VARIABLES FOR MACROSCOPIC FIELDS AND DISTRIBUTION FUNCTIONS             ####
###############################################################################################################

# initializing the electric and magnetic fields
def initialize_field(Ny=10, Nx=10):
    return np.zeros((Ny, Nx), dtype=np.float32, order='C')

Ex_inc, Ey_inc, Ez_inc, Hx_inc, Hy_inc, Hz_inc = [initialize_field(Ny, Nx) for _ in range(6)]
Ex_tot, Ey_tot, Ez_tot, Hx_tot, Hy_tot, Hz_tot = [initialize_field(Ny, Nx) for _ in range(6)]
Ex_scat, Ey_scat, Ez_scat, Hx_scat, Hy_scat, Hz_scat = [initialize_field(Ny, Nx) for _ in range(6)]


# initializing the distribution functions of electric and magnetic fields
def initilize_dis_func(Ny=10, Nx=10, Q=7):
    return np.zeros((Ny, Nx, Q), dtype=np.float32, order='C')

ex_inc, ey_inc, ez_inc, hx_inc, hy_inc, hz_inc = [initilize_dis_func(Ny, Nx, Q) for _ in range(6)]
ex_tot, ey_tot, ez_tot, hx_tot, hy_tot, hz_tot = [initilize_dis_func(Ny, Nx, Q) for _ in range(6)]


# initilizing the domain properties
def initialize_material_properties(er1=1, mur1=1, Ny=10, Nx=10):
    yield np.ones((Ny, Nx), dtype=np.float32, order='C') * er1
    yield np.ones((Ny, Nx), dtype=np.float32, order='C') * mur1
    return

er_inc, mur_inc = initialize_material_properties(er1, mur1, Ny, Nx)
er_tot, mur_tot = initialize_material_properties(er1, mur1, Ny, Nx)
###############################################################################################################




################################################################
###                        SCATTERER                        ####
################################################################

a1 = a

# initilizing the polar coordinates
r, phi = [initialize_field(Ny, Nx) for _ in range(2)]

# center of the scatterer
cx = Nx//2 + 0.5
cy = Ny//2 + 0.5

# converting to polar coordinates
Module_Geometry.carToPolar(r, phi, Ny, Nx, cy, cx)

# scatterer particle
scatterer = Module_Geometry.circle(r, a, Ny, Nx)

er_tot[scatterer] = er2
################################################################







t1 = time.time()



fig = plt.figure(figsize = (3.35*3, 3.35), dpi=600)
gs=GridSpec(1,2)
plt.ion()

fps = 10

Video_Name = "EzTotScat.mp4"

FFMpegWriter = animation.writers['ffmpeg']

metadata = dict(title='Ez_tot', artist='', comment='Movie support!')
writer = FFMpegWriter(fps=fps, metadata=metadata)

with writer.saving(fig, Video_Name,300):
    for t in range(int(Time)):
        
        

        #################################################################################################################
        ########                                         LBM CALCULATION                                          #######
        #################################################################################################################

        # initialization of macroscopic fields
        myclib.initializeField(Ex_inc, Ey_inc, Ez_inc, Hx_inc, Hy_inc, Hz_inc, Ny, Nx)
        myclib.initializeField(Ex_tot, Ey_tot, Ez_tot, Hx_tot, Hy_tot, Hz_tot, Ny, Nx)

        # computation of macroscopic fields from distribution function
        myclib.macroField(ex_inc, er_inc, Ex_inc, Ny, Nx, Q)
        myclib.macroField(ey_inc, er_inc, Ey_inc, Ny, Nx, Q)
        myclib.macroField(ez_inc, er_inc, Ez_inc, Ny, Nx, Q)

        myclib.macroField(hx_inc, mur_inc, Hx_inc, Ny, Nx, Q)
        myclib.macroField(hy_inc, mur_inc, Hy_inc, Ny, Nx, Q)
        myclib.macroField(hz_inc, mur_inc, Hz_inc, Ny, Nx, Q)


        myclib.macroField(ex_tot, er_tot, Ex_tot, Ny, Nx, Q)
        myclib.macroField(ey_tot, er_tot, Ey_tot, Ny, Nx, Q)
        myclib.macroField(ez_tot, er_tot, Ez_tot, Ny, Nx, Q)

        myclib.macroField(hx_tot, mur_tot, Hx_tot, Ny, Nx, Q)
        myclib.macroField(hy_tot, mur_tot, Hy_tot, Ny, Nx, Q)
        myclib.macroField(hz_tot, mur_tot, Hz_tot, Ny, Nx, Q)


        

        
        if (t >= 0):
            
            # source wave
            planeWaveTM(Ez_inc, Hy_inc, t, omega, xloc, ymin, ymax)
            planeWaveTM(Ez_tot, Hy_tot, t, omega, xloc, ymin, ymax)

            # calculation of scattered fields
            Ez_scat = Ez_tot - Ez_inc
            Hx_scat = Hx_tot - Hx_inc
            Hy_scat = Hy_tot - Hy_inc

            Ez_scat[scatterer] = 0
            Hx_scat[scatterer] = 0
            Hy_scat[scatterer] = 0


            # collision and streaming (the 2 steps of LBM) when field is forced
            myclib.collStreamForcingNode(ex_inc, ey_inc, ez_inc, hx_inc, hy_inc, hz_inc, Ex_inc, Ey_inc, Ez_inc, Hx_inc, Hy_inc, Hz_inc, er_inc, mur_inc, Ny, Nx, Q, xloc, ymin, ymax)
            myclib.collStreamForcingNode(ex_tot, ey_tot, ez_tot, hx_tot, hy_tot, hz_tot, Ex_tot, Ey_tot, Ez_tot, Hx_tot, Hy_tot, Hz_tot, er_tot, mur_tot, Ny, Nx, Q, xloc, ymin, ymax)
            
        else:

            # calculation of scattered fields
            Ez_scat = Ez_tot - Ez_inc
            Hx_scat = Hx_tot - Hx_inc
            Hy_scat = Hy_tot - Hy_inc

            Ez_scat[scatterer] = 0
            Hx_scat[scatterer] = 0
            Hy_scat[scatterer] = 0

            
            # collision and streaming (the 2 steps of LBM) when field is not forced
            myclib.collStream(ex_inc, ey_inc, ez_inc, hx_inc, hy_inc, hz_inc, Ex_inc, Ey_inc, Ez_inc, Hx_inc, Hy_inc, Hz_inc, er_inc, mur_inc, Ny, Nx, Q)
            myclib.collStream(ex_tot, ey_tot, ez_tot, hx_tot, hy_tot, hz_tot, Ex_tot, Ey_tot, Ez_tot, Hx_tot, Hy_tot, Hz_tot, er_tot, mur_tot, Ny, Nx, Q)
            
        ###############################################################################################################




        

    ###############################################################################################################
    ##########                                           ANIMATION                                       ##########
    ###############################################################################################################

        if (t%10 == 0):
                
            fig.clear()

            ########### Plot parameters ################
            plt.rc('font', family = 'serif', size = 10)
            plt.rc('xtick', labelsize = 10)
            plt.rc('ytick', labelsize = 10)
            plt.rc('lines', markersize = 2, lw = 0.75)
            plt.rc('text', usetex = True)
            plt.rcParams['image.cmap']='seismic'
            ###########################################
            
            
                        
            ax1 = fig.add_subplot(gs[0,0])
            plt.title(r'$E_z^{tot}$')
            im1 = plt.imshow(Ez_tot, vmin = -1, vmax = 1, cmap='seismic', origin='lower')
            im = plt.imshow(er_tot, extent=(0, Nx, 0, Ny), cmap='binary', origin='lower', alpha=0.1)
            
            ax1.set_xticks(np.linspace(0,Nx,4))
            ax1.set_xticklabels([])
            ax1.set_xlabel(r'$x$', fontsize=12)
            ax1.set_yticks(np.linspace(0,Ny,4))
            ax1.set_yticklabels([])
            ax1.set_ylabel(r'$y$', fontsize=12)


        
            ax2 = fig.add_subplot(gs[0,1])
            plt.title(r'$E_z^{scat}$')           
            im2 = plt.imshow(Ez_scat, vmin = -1, vmax = 1, cmap='seismic', origin='lower')
            plt.imshow(er_tot, extent=(0, Nx, 0, Ny), cmap='binary', origin='lower', alpha=0.1)
            
            ax2.set_xticks(np.linspace(0,Nx,4))
            ax2.set_xticklabels([])
            ax2.set_xlabel(r'$x$', fontsize=12)
            ax2.set_yticks(np.linspace(0,Ny,4))
            ax2.set_yticklabels([])
##            ax2.set_ylabel(r'$y$', fontsize=12)

            
       
##            plt.savefig(pictures+"/pic."+str(t).zfill(4)+".png")
            writer.grab_frame()     

    ###############################################################################################################
    ###############################################################################################################



###############################################################################################################
        t2 = time.time()
        
        if (t > 0 and t%100 == 0):
            remaining_time = (t2 - t1) * (int(Time) - t) / (t*60)
            print(f"Approximate time left: {remaining_time:.2f} minutes", end="\r")

t3 = time.time()
total_time = (t3 - t0) / 60
print(f"\nTotal time taken: {total_time:.2f} minutes")
###############################################################################################################
