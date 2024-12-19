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


from ModuleSourceOfEmWave import planeWaveTM
import ModuleGeometry


t0 = time.time()


pictures = 'gallery'
if not os.path.exists(pictures):
    os.makedirs(pictures)


directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)




# Accessing command line arguments
parameters = sys.argv


######################


a, ratio = 40, 1   # ratio = a / wavelength
n = 20
Nx, Ny, Nz = n*a, n*a, 1  # size of the computational domain

er1, mur1, er2 = 1, 1, 10000   # material properties i.e. permittivity and permeabilty


#############################
# boundary of EM wave source
xloc = 0
ymin = 0
ymax = Ny
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



noOfPeriods = 1000
noOfReflections = 0

Time = 3 * (Nx * np.sqrt(er1) + noOfReflections * 2 * a * np.sqrt(er2)) + noOfPeriods * period

print("Number of time steps :", int(Time))




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

scatterer = []


################################################################
###                        SCATTERER 1                      ####
################################################################

a1 = a


for i in range(Ny):
    for j in range(Nx):
        for n in range(0, 10, 2):
            for m in range(0, 10, 2):
                cy = 2*a*n + 2*a -1
                cx = 2*a*m + 2*a -1

##                print(cx, cy)

                Y = i - cy
                X = j - cx
                if (X <= a and X >= -a and Y <= a and Y >= -a):
                    er[0, i, j] = er2
            
        
##
##        # scatterer particle
##        scatterer.append(ModuleGeometry.square(a, Nz, Ny, Nx, cy, cx))
##        print(scatterer)
##
##er[scatterer] = er2
################################################################





###############################################################################################################
########                                           SHARED LIBRARY                                   ###########
###############################################################################################################

# loading the shared file (c library)
path = os.getcwd()
myclib = CDLL(os.path.join(path, "LBM.so"))

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


myclib.collForcingNode.argtypes = [P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P3D, P3D, P3D, P3D, P3D, P3D, P3D, P3D, c_int, c_int, c_int, c_int, c_int, c_int, c_int]
myclib.collForcingNode.restype  = None


myclib.streaming.argtypes = [P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, P4D, c_int, c_int, c_int, c_int]
myclib.streaming.restype  = None

###############################################################################################################


U = np.zeros(int(Time))


t1 = time.time()



fig = plt.figure(figsize = (3.35*1, 2.85), dpi=600)
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
            planeWaveTM(Ez, Hy, t, omega, xloc, ymin, ymax)

            # collision and streaming (the 2 steps of LBM) when field is forced
            myclib.collForcingNode(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Nz, Ny, Nx, Q, xloc, ymin, ymax)
            myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Nz, Ny, Nx, Q)


        else:
            
            # collision and streaming (the 2 steps of LBM) when field is not forced
            myclib.collNotForcingNode(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Nz, Ny, Nx, Q)
            myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Nz, Ny, Nx, Q)
            
        ###############################################################################################################


        U[t] = np.sum(0.5 * (er1 * (Ex**2 + Ey**2 + Ez**2) + mur1 * (Hx**2 + Hy**2 + Hz**2)))

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
            
            nx = Nx//(4*a)
            ny = Ny//(4*a)
                        
            ax1 = fig.add_subplot(gs[0,0])
            plt.title(r'$E_z^{tot}$')
            im1 = plt.imshow(Ez[0], vmin = -1, vmax = 1, cmap='seismic', origin='lower')
##            im = plt.imshow(er[0], extent=(0, Nx, 0, Ny), cmap='binary', origin='lower', alpha=0.5)


            # boundary of the hexagon
##            n = 50
##            plt.plot(np.ones(n)*(cx-a), np.linspace(cx-a, cx+a, n), 'k-')
##            plt.plot(np.ones(n)*(cx+a), np.linspace(cx-a, cx+a, n), 'k-')
##            plt.plot(np.linspace(cx-a, cx, n), np.linspace(cx+a, cx+2*a, n), 'k-')
##            plt.plot(np.linspace(cx, cx+a, n), np.linspace(cx+2*a, cx+a, n), 'k-')
##            plt.plot(np.linspace(cx-a, cx, n), np.linspace(cx-a, cx-2*a, n), 'k-')
##            plt.plot(np.linspace(cx, cx+a, n), np.linspace(cx-2*a, cx-a, n), 'k-')

            
            ax1.set_xticks(np.linspace(0,Nx,5))
            ax1.set_xticklabels([0, nx, 2*nx, 3*nx, 4*nx])
            ax1.set_xlabel(r'$x/a$', fontsize=12)
            ax1.set_yticks(np.linspace(0,Ny,5))
            ax1.set_yticklabels([0, ny, 2*ny, 3*ny, 4*ny])
            ax1.set_ylabel(r'$y/a$', fontsize=12)

            fig.colorbar(im1, ax = ax1, location='right', shrink=0.99, aspect=25)

            plt.tight_layout()
       
            plt.savefig(pictures+"/pic."+str(t).zfill(4)+".svg")
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






energy = open(directory+"/energy_ratio_{}_a_{}.txt".format(ratio, a), "w")
np.savetxt(energy, U, fmt='%.4e')
energy.close()
