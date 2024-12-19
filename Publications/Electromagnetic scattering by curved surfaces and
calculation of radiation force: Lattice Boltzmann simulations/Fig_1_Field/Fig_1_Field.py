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




# Accessing command line arguments
parameters = sys.argv


######################


a, ratio = 25, 1.0   # ratio = a / wavelength
n = 20
Nx, Ny, Nz = n*a, n*a, 1  # size of the computational domain

er1, mur1, er2 = 1, 1, 10   # material properties i.e. permittivity and permeabilty


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



noOfPeriods = 0
noOfReflections = 0

Time = 3 * (Nx * np.sqrt(er1) + noOfReflections * 2 * a * np.sqrt(er2)) + noOfPeriods * period

print("Number of time steps :", int(Time))




###############################################################################################################
####      DEFINING AND INITILIZING VARIABLES FOR MACROSCOPIC FIELDS AND DISTRIBUTION FUNCTIONS             ####
###############################################################################################################

# initializing the electric and magnetic fields
def initialize_field(Nz=1, Ny=10, Nx=10):
    return np.zeros((Nz, Ny, Nx), dtype=np.float32, order='C')

Ex_inc, Ey_inc, Ez_inc, Hx_inc, Hy_inc, Hz_inc = [initialize_field(Nz, Ny, Nx) for _ in range(6)]
Ex_tot, Ey_tot, Ez_tot, Hx_tot, Hy_tot, Hz_tot = [initialize_field(Nz, Ny, Nx) for _ in range(6)]
Ex_scat, Ey_scat, Ez_scat, Hx_scat, Hy_scat, Hz_scat = [initialize_field(Nz, Ny, Nx) for _ in range(6)]


# initializing the distribution functions of electric and magnetic fields
def initilize_dis_func(Nz=1, Ny=10, Nx=10, Q=7):
    return np.zeros((Nz, Ny, Nx, Q), dtype=np.float32, order='C')

ex_inc, ey_inc, ez_inc, hx_inc, hy_inc, hz_inc = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]
exb_inc, eyb_inc, ezb_inc, hxb_inc, hyb_inc, hzb_inc = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]
ex_tot, ey_tot, ez_tot, hx_tot, hy_tot, hz_tot = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]
exb_tot, eyb_tot, ezb_tot, hxb_tot, hyb_tot, hzb_tot = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]


# initilizing the domain properties
def initialize_material_properties(er1=1, mur1=1, Nz=1, Ny=10, Nx=10):
    yield np.ones((Nz, Ny, Nx), dtype=np.float32, order='C') * er1
    yield np.ones((Nz, Ny, Nx), dtype=np.float32, order='C') * mur1
    return

er_inc, mur_inc = initialize_material_properties(er1, mur1, Nz, Ny, Nx)
er_tot, mur_tot = initialize_material_properties(er1, mur1, Nz, Ny, Nx)
###############################################################################################################




################################################################
###                        SCATTERER 1                      ####
################################################################

a = a

# initilizing the polar coordinates
r, phi = [initialize_field(Nz, Ny, Nx) for _ in range(2)]

# center of the scatterer
cx = Nx//2 + 0.5
cy = Ny//2 + 0.5

# converting to polar coordinates
ModuleGeometry.carToPolar(r, phi, Nz, Ny, Nx, cy, cx)

# scatterer particle
scatterer = ModuleGeometry.circle(r, a, Nz, Ny, Nx)

er_tot[scatterer] = er2
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


t1 = time.time()




fig = plt.figure(figsize = (6.69, 2.75), dpi=600, tight_layout=True)
gs=GridSpec(1,3)
plt.ion()



fps = 20

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
        myclib.initializeField(Ex_inc, Ey_inc, Ez_inc, Hx_inc, Hy_inc, Hz_inc, Nz, Ny, Nx)
        myclib.initializeField(Ex_tot, Ey_tot, Ez_tot, Hx_tot, Hy_tot, Hz_tot, Nz, Ny, Nx)

        # computation of macroscopic fields from distribution function
        myclib.macroField(ex_inc, er_inc, Ex_inc, Nz, Ny, Nx, Q)
        myclib.macroField(ey_inc, er_inc, Ey_inc, Nz, Ny, Nx, Q)
        myclib.macroField(ez_inc, er_inc, Ez_inc, Nz, Ny, Nx, Q)

        myclib.macroField(hx_inc, mur_inc, Hx_inc, Nz, Ny, Nx, Q)
        myclib.macroField(hy_inc, mur_inc, Hy_inc, Nz, Ny, Nx, Q)
        myclib.macroField(hz_inc, mur_inc, Hz_inc, Nz, Ny, Nx, Q)


        myclib.macroField(ex_tot, er_tot, Ex_tot, Nz, Ny, Nx, Q)
        myclib.macroField(ey_tot, er_tot, Ey_tot, Nz, Ny, Nx, Q)
        myclib.macroField(ez_tot, er_tot, Ez_tot, Nz, Ny, Nx, Q)

        myclib.macroField(hx_tot, mur_tot, Hx_tot, Nz, Ny, Nx, Q)
        myclib.macroField(hy_tot, mur_tot, Hy_tot, Nz, Ny, Nx, Q)
        myclib.macroField(hz_tot, mur_tot, Hz_tot, Nz, Ny, Nx, Q)


        

        # plane wave enforcement at the left boundary of the domain
        if (t >= 0):
            
            # source wave
            planeWaveTM(Ez_inc, Hy_inc, t, omega, xloc, ymin, ymax)
            planeWaveTM(Ez_tot, Hy_tot, t, omega, xloc, ymin, ymax)

            # calculation of scattered fields
            Ez_scat = Ez_tot - Ez_inc
            Hx_scat = Hx_tot - Hx_inc
            Hy_scat = Hy_tot - Hy_inc

            Ez_scat[scatterer] = 0
            

            # collision and streaming (the 2 steps of LBM) when field is forced
            myclib.collForcingNode(ex_inc, ey_inc, ez_inc, hx_inc, hy_inc, hz_inc, exb_inc, eyb_inc, ezb_inc, hxb_inc, hyb_inc, hzb_inc, Ex_inc, Ey_inc, Ez_inc, Hx_inc, Hy_inc, Hz_inc, er_inc, mur_inc, Nz, Ny, Nx, Q, xloc, ymin, ymax)
            myclib.collForcingNode(ex_tot, ey_tot, ez_tot, hx_tot, hy_tot, hz_tot, exb_tot, eyb_tot, ezb_tot, hxb_tot, hyb_tot, hzb_tot, Ex_tot, Ey_tot, Ez_tot, Hx_tot, Hy_tot, Hz_tot, er_tot, mur_tot, Nz, Ny, Nx, Q, xloc, ymin, ymax)

            myclib.streaming(ex_inc, ey_inc, ez_inc, hx_inc, hy_inc, hz_inc, exb_inc, eyb_inc, ezb_inc, hxb_inc, hyb_inc, hzb_inc, Nz, Ny, Nx, Q)
            myclib.streaming(ex_tot, ey_tot, ez_tot, hx_tot, hy_tot, hz_tot, exb_tot, eyb_tot, ezb_tot, hxb_tot, hyb_tot, hzb_tot, Nz, Ny, Nx, Q)
        else:

            # calculation of scattered fields
            Ez_scat = Ez_tot - Ez_inc
            Hx_scat = Hx_tot - Hx_inc
            Hy_scat = Hy_tot - Hy_inc

            Ez_scat[scatterer] = 0
            

            
            # collision and streaming (the 2 steps of LBM) when field is not forced
            myclib.collStream(ex_inc, ey_inc, ez_inc, hx_inc, hy_inc, hz_inc, exb_inc, eyb_inc, ezb_inc, hxb_inc, hyb_inc, hzb_inc, Ex_inc, Ey_inc, Ez_inc, Hx_inc, Hy_inc, Hz_inc, er_inc, mur_inc, Nz, Ny, Nx, Q)
            myclib.collStream(ex_tot, ey_tot, ez_tot, hx_tot, hy_tot, hz_tot, exb_tot, eyb_tot, ezb_tot, hxb_tot, hyb_tot, hzb_tot, Ex_tot, Ey_tot, Ez_tot, Hx_tot, Hy_tot, Hz_tot, er_tot, mur_tot, Nz, Ny, Nx, Q)

            myclib.streaming(ex_inc, ey_inc, ez_inc, hx_inc, hy_inc, hz_inc, exb_inc, eyb_inc, ezb_inc, hxb_inc, hyb_inc, hzb_inc, Nz, Ny, Nx, Q)
            myclib.streaming(ex_tot, ey_tot, ez_tot, hx_tot, hy_tot, hz_tot, exb_tot, eyb_tot, ezb_tot, hxb_tot, hyb_tot, hzb_tot, Nz, Ny, Nx, Q)
        ###############################################################################################################


        
        

    ###############################################################################################################
    ##########                                           ANIMATION                                       ##########
    ###############################################################################################################

        if (t%5 == 0):
                
            fig.clear()

            ########### Plot parameters ################
            plt.rc('font', family = 'serif', size = 8)
            plt.rc('xtick', labelsize = 8)
            plt.rc('ytick', labelsize = 8)
            plt.rc('lines', markersize = 3, lw = 0.5)
            plt.rc('text', usetex = True)
            plt.rcParams['image.cmap']='seismic'
            ###########################################

            nx = Nx//(5*a)
            ny = Ny//(5*a)

            # incident field
            Ez_inc[0, :, 200:] = 0
            
            ax0 = fig.add_subplot(gs[0,0])
            plt.title(r'Incident field ($\mathcal{E}_z^I$)')
            plt.imshow(Ez_inc[0], vmin = -1, vmax = 1, cmap='seismic', origin='lower')
            plt.imshow(er_tot[0], extent=(0, Nx, 0, Ny), cmap='binary', origin='lower', alpha=0.3)
            
            ax0.set_xticks(np.linspace(0,Nx,6))
            ax0.set_xticklabels([0, nx, 2*nx, 3*nx, 4*nx, 5*nx])
            ax0.set_xlabel(r'$x/a$', fontsize=10)
            ax0.set_yticks(np.linspace(0,Ny,6))
            ax0.set_yticklabels([0, ny, 2*ny, 3*ny, 4*ny, 5*ny])
            ax0.set_ylabel(r'$y/a$', fontsize=10)

            ax0.plot([268,250], [268,250], lw=0.5, color='k')
                
            ax0.arrow(232+63*1.58,269, -60, 0, lw=0.5, head_width=15, length_includes_head=True, edgecolor='k', facecolor='k')
            plt.text(232+65,268+25, r'$a$', color='k', bbox={'fc':'white', 'ec':'white','boxstyle':'square, pad=0.1'})
                
            ax0.annotate('', xy=(velocity*t - 300 - wavelength, 400), xytext=(velocity*t - 300 - 4*wavelength, 400), arrowprops=dict(fc="white", ec="white", width=0.1, headwidth=3, headlength=5, shrink=0.05))
            ax0.annotate('', xy=(velocity*t - 300, 400), xytext=(velocity*t - 300 + 2*wavelength, 400), arrowprops=dict(fc="k", ec="k", width=0.1, headwidth=3, headlength=5, shrink=0.05))
            plt.text(velocity*t - 300  - 3.75*wavelength, 420, r'$\lambda$', color='k', bbox={'fc':'white', 'ec':'white','boxstyle':'square, pad=0.1'})
                
            ax0.annotate('', xy=(200 + 100, 150), xytext=(200, 150), arrowprops=dict(fc="k", ec="k", width=3, headwidth=6, headlength=8, shrink=0.05))
            ax0.text(200 + 30, 165, r'$v$')
            ax0.text(0.05*Nx, 0.9*Ny, r'$(a)$', color='k', bbox={'fc':'white', 'ec':'white','boxstyle':'square, pad=0.1'})


            
            # total field        
            ax1 = fig.add_subplot(gs[0,1])
            plt.title(r'Total field ($\mathcal{E}_z^{total}$)')
            im1 = plt.imshow(Ez_tot[0], vmin = -1, vmax = 1, cmap='seismic', origin='lower')
            im = plt.imshow(er_tot[0], extent=(0, Nx, 0, Ny), cmap='binary', origin='lower', alpha=0.3)
            
            ax1.set_xticks(np.linspace(0,Nx,6))
            ax1.set_xticklabels([0, nx, 2*nx, 3*nx, 4*nx, 5*nx])
            ax1.set_xlabel(r'$x/a$', fontsize=10)
            ax1.set_yticks(np.linspace(0,Ny,6))
            ax1.set_yticklabels([])

            ax1.text(0.05*Nx, 0.9*Ny, r'$(b)$', color='k', bbox={'fc':'white', 'ec':'white','boxstyle':'square, pad=0.1'})


            # scattered field
            ax2 = fig.add_subplot(gs[0,2])
            plt.title(r'Scattered field ($\mathcal{E}_z^S$)')
            im2 = plt.imshow(Ez_scat[0], vmin = -1, vmax = 1, cmap='seismic', origin='lower')
            plt.imshow(er_tot[0], extent=(0, Nx, 0, Ny), cmap='binary', origin='lower', alpha=0.3)
            
            ax2.set_xticks(np.linspace(0,Nx,6))
            ax2.set_xticklabels([0, nx, 2*nx, 3*nx, 4*nx, 5*nx])
            ax2.set_xlabel(r'$x/a$', fontsize=10)
            ax2.set_yticks(np.linspace(0,Ny,6))
            ax2.set_yticklabels([])

            ax2.plot([100,400], [250,250], lw=0.5, color='k', linestyle='--')
            ax2.annotate('', xy=(318, 368), xytext=(250, 250), arrowprops=dict(fc="k", ec="k", width=0.1, headwidth=3, headlength=5))
            ax2.plot(325,375, 'ok', markersize=1)
            ax2.text(305, 275, r'$\mathbf{\phi}$', bbox={'fc':'white', 'ec':'white','boxstyle':'square, pad=0.1'})
            ax2.text(255, 340, r'$r$', color='k', bbox={'fc':'white', 'ec':'white','boxstyle':'square, pad=0.1'})
            ax2.text(345, 370, r'$P$', color='k', bbox={'fc':'white', 'ec':'white','boxstyle':'square, pad=0.1'})
            ax2.text(0.05*Nx, 0.9*Ny, r'$(c)$', color='k', bbox={'fc':'white', 'ec':'white','boxstyle':'square, pad=0.1'})


            
       
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
