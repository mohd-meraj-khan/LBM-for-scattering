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

from ModuleTraction import tractionFx, tractionFy, tractionFz
from ModuleSourceOfEmWave import planeWaveTM
import ModuleGeometry




pictures = 'gallery'
if not os.path.exists(pictures):
    os.makedirs(pictures)




# Accessing command line arguments
parameters = sys.argv


######################


a, ratio = 25, 0.5   # ratio = a / wavelength
n = 4
Nx, Ny, Nz = n*a, n*a, 1  # size of the computational domain

er1, mur1, er2, er3 = 1, 1, 2, 5   # material properties i.e. permittivity and permeabilty


#############################
# boundary of EM wave source
xmin = 0
xmax = 0 + 1
ymin = 0
ymax = Ny + 1 
#############################


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



noOfPeriods = 1
noOfReflections = 10

Time = 3 * (Nx * np.sqrt(er1) + noOfReflections * 2 * a * np.sqrt(er2)) + noOfPeriods * period

print(Time)




###############################################################################################################
####      DEFINING AND INITILIZING VARIABLES FOR MACROSCOPIC FIELDS AND DISTRIBUTION FUNCTIONS             ####
###############################################################################################################

# initializing the electric and magnetic fields
def initialize_field(Nz=1, Ny=10, Nx=10):
    return np.zeros((Nz, Ny, Nx), dtype=np.float32, order='C')

Ex, Ey, Ez, Hx, Hy, Hz = [initialize_field(Nz, Ny, Nx) for _ in range(6)]


# initializing the distribution functions of electric and magnetic fields
def initilize_dis_func(Nz=1, Ny=10, Nx=10, Q=7):
    return np.zeros((Nz, Ny, Nx, Q), dtype=np.float32, order='C')

ex, ey, ez, hx, hy, hz = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]
exb, eyb, ezb, hxb, hyb, hzb = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]


# initilizing the domain properties
def initialize_material_properties(er1=1, mur1=1, Nz=1, Ny=10, Nx=10):
    yield np.ones((Nz, Ny, Nx), dtype=np.float32, order='C') * er1
    yield np.ones((Nz, Ny, Nx), dtype=np.float32, order='C') * mur1
    return

er, mur = initialize_material_properties(er1, mur1, Nz, Ny, Nx)
###############################################################################################################




################################################################
###                        SCATTERER 1                      ####
################################################################

a1 = a

# initilizing the polar coordinates
r1, phi1 = [initialize_field(Nz, Ny, Nx) for _ in range(2)]

# center of the scatterer
cx1 = Nx//2 + 0.5
cy1 = Ny//2 + 0.5

# converting to polar coordinates
ModuleGeometry.carToPolar(r1, phi1, Nz, Ny, Nx, cy1, cx1)

# scatterer particle
ModuleGeometry.circle(r1, phi1, er, er2, a1, Nz, Ny, Nx)
################################################################





###############################################################################################################
########                                           SHARED LIBRARY                                   ###########
###############################################################################################################

# loading the shared file (c library)
path = os.getcwd()
myclib = CDLL(os.path.join(path, "LBM_OpenMP.so"))

# defining 3D and 4D pointers (LBM runs in C, for that pointer is needed)
P3D = np.ctypeslib.ndpointer(dtype=np.float32, ndim=3, flags="C")
P4D = np.ctypeslib.ndpointer(dtype=np.float32, ndim=4, flags="C")

# calculation of macroscopic fields (FUNCTION PROTOTYPE)
myclib.macroField.argtypes = [P4D, P3D, P3D, c_int, c_int, c_int, c_int]
myclib.macroField.restype  = None

# initilization of macroscopic fields (FUNCTION PROTOTYPE)
myclib.initializeField.argtypes = [P3D, P3D, P3D, P3D, P3D, P3D, c_int, c_int, c_int]
myclib.initializeField.restype  = None

# collision + streaming (FUNCTION PROTOTYPE)
myclib.collNotForcingNode.argtypes = [P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P3D, P3D, P3D, P3D, P3D, P3D, P3D, P3D, c_int, c_int, c_int, c_int]
myclib.collNotForcingNode.restype  = None


myclib.collForcingNode.argtypes = [P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P3D, P3D, P3D, P3D, P3D, P3D, P3D, P3D, c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int]
myclib.collForcingNode.restype  = None


myclib.streaming.argtypes = [P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, c_int, c_int, c_int, c_int]
myclib.streaming.restype  = None

###############################################################################################################





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
        
        t1 = time.time()

        #################################################################################################################
        ########                                         LBM CALCULATION                                          #######
        #################################################################################################################

        # initialization of macroscopic fields
        myclib.initializeField(Ex, Ey, Ez, Hx, Hy, Hz, Nz, Ny, Nx)

        # computation of macroscopic fields from distribution function
        myclib.macroField(ex, er, Ex, Nz, Ny, Nx, Q)
        myclib.macroField(ey, er, Ey, Nz, Ny, Nx, Q)
        myclib.macroField(ez, er, Ez, Nz, Ny, Nx, Q)

        myclib.macroField(hx, mur, Hx, Nz, Ny, Nx, Q)
        myclib.macroField(hy, mur, Hy, Nz, Ny, Nx, Q)
        myclib.macroField(hz, mur, Hz, Nz, Ny, Nx, Q)

        
        if (t >= 0):
            
            # source wave
            planeWaveTM(Ez, Hy, t, omega, xmin, xmax, ymin, ymax)

            # collision and streaming (the 2 steps of LBM) when field is forced
            myclib.collForcingNode(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Nz, Ny, Nx, Q, xmin, xmax, ymin, ymax)
            myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Nz, Ny, Nx, Q)


        else:
            
            # collision and streaming (the 2 steps of LBM) when field is not forced
            myclib.collNotForcingNode(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Nz, Ny, Nx, Q)
            myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Nz, Ny, Nx, Q)
            
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
            im1 = plt.imshow(Ez[0], vmin = -1, vmax = 1, cmap='seismic', origin='lower')
            im = plt.imshow(er[0], extent=(0, Nx, 0, Ny), cmap='binary', origin='lower', alpha=0.1)
            
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



