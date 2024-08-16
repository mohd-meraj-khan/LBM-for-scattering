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

Ex, Ey, Ez, Hx, Hy, Hz = [initialize_field(Ny, Nx) for _ in range(6)]


# initializing the distribution functions of electric and magnetic fields
def initilize_dis_func(Ny=10, Nx=10, Q=7):
    return np.zeros((Ny, Nx, Q), dtype=np.float32, order='C')

ex, ey, ez, hx, hy, hz = [initilize_dis_func(Ny, Nx, Q) for _ in range(6)]
exb, eyb, ezb, hxb, hyb, hzb = [initilize_dis_func(Ny, Nx, Q) for _ in range(6)]


# initilizing the domain properties
def initialize_material_properties(er1=1, mur1=1, Ny=10, Nx=10):
    yield np.ones((Ny, Nx), dtype=np.float32, order='C') * er1
    yield np.ones((Ny, Nx), dtype=np.float32, order='C') * mur1
    return

er, mur = initialize_material_properties(er1, mur1, Ny, Nx)
###############################################################################################################




################################################################
###                        SCATTERER                        ####
################################################################


# initilizing the polar coordinates
r, phi = [initialize_field(Ny, Nx) for _ in range(2)]

# center of the scatterer
cx = Nx//2 + 0.5
cy = Ny//2 + 0.5

# converting to polar coordinates
Module_Geometry.carToPolar(r, phi, Ny, Nx, cy, cx)

# scatterer particle
scatterer = Module_Geometry.circle(r, a, Ny, Nx)

er[scatterer] = er2
################################################################






t1 = time.time()


fig = plt.figure(figsize = (3.35*1, 3.35), dpi=600)
gs=GridSpec(1,1)
plt.ion()

fps = 10

Video_Name = "Ez_tot.mp4"

FFMpegWriter = animation.writers['ffmpeg']

metadata = dict(title='Ez_tot', artist='', comment='Movie support!')
writer = FFMpegWriter(fps=fps, metadata=metadata)

with writer.saving(fig, Video_Name,300):
    for t in range(int(Time)):
        
        

        #################################################################################################################
        ########                                         LBM CALCULATION                                          #######
        #################################################################################################################

        # initialization of macroscopic fields
        myclib.initializeField(Ex, Ey, Ez, Hx, Hy, Hz, Ny, Nx)

        # computation of macroscopic fields from distribution function
        myclib.macroField(ex, er, Ex, Ny, Nx, Q)
        myclib.macroField(ey, er, Ey, Ny, Nx, Q)
        myclib.macroField(ez, er, Ez, Ny, Nx, Q)

        myclib.macroField(hx, mur, Hx, Ny, Nx, Q)
        myclib.macroField(hy, mur, Hy, Ny, Nx, Q)
        myclib.macroField(hz, mur, Hz, Ny, Nx, Q)

        
        if (t >= 0):
            
            # source wave
            planeWaveTM(Ez, Hy, t, omega, xloc, ymin, ymax)

            # collision and streaming (the 2 steps of LBM) when field is forced
            myclib.collForcingNode(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Ny, Nx, Q, xloc, ymin, ymax)
            myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ny, Nx, Q)


        else:
            
            # collision and streaming (the 2 steps of LBM) when field is not forced
            myclib.collNotForcingNode(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Ny, Nx, Q)
            myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ny, Nx, Q)
            
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
            im1 = plt.imshow(Ez, vmin = -1, vmax = 1, cmap='seismic', origin='lower')
            im = plt.imshow(er, extent=(0, Nx, 0, Ny), cmap='binary', origin='lower', alpha=0.1)
            
            ax1.set_xticks(np.linspace(0,Nx,4))
            ax1.set_xticklabels([])
            ax1.set_xlabel(r'$x$', fontsize=12)
            ax1.set_yticks(np.linspace(0,Ny,4))
            ax1.set_yticklabels([])
            ax1.set_ylabel(r'$y$', fontsize=12)


       
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
