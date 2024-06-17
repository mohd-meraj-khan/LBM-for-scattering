from ctypes import CDLL, POINTER
from ctypes import c_size_t, c_double, c_float, c_int
import numpy as np
from scipy.interpolate import RectBivariateSpline
import sys
import os
import time

from ModuleTraction import tractionFx, tractionFy, tractionFz
from ModuleSourceOfEmWave import planeWaveTM
import ModuleGeometry




directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)




# Accessing command line arguments
parameters = sys.argv


######################


a, ratio = 25, 0.95   # ratio = a / wavelength
n = 4
Nx, Ny, Nz = n*a, n*a, 1  # size of the computational domain

er1, mur1, er2, er3 = 1, 1, 4, 5   # material properties i.e. permittivity and permeabilty


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
noOfReflections = 30

Time = 3 * (Nx * np.sqrt(er1) + noOfReflections * 2 * a * np.sqrt(er2)) + noOfPeriods * period

print(Time)




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
ex_tot, ey_tot, ez_tot, hx_tot, hy_tot, hz_tot = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]


# initilizing the domain properties
def initialize_material_properties(er1=1, mur1=1, Nz=1, Ny=10, Nx=10):
    yield np.ones((Nz, Ny, Nx), dtype=np.float32, order='C') * er1
    yield np.ones((Nz, Ny, Nx), dtype=np.float32, order='C') * mur1
    return

er_inc, mur_inc = initialize_material_properties(er1, mur1, Nz, Ny, Nx)
er_tot, mur_tot = initialize_material_properties(er1, mur1, Nz, Ny, Nx)
###############################################################################################################




################################################################
################### for interpolation ##########################

dtheta = 0.1
theta = np.arange(0, 360 + 0.01, dtheta)

# surface normal (the integration is done at a circle enclosing the scatterer)
nx = np.cos(theta*np.pi/180)
ny = np.sin(theta*np.pi/180)
nz = 0

# coordinates of the all grids in the domain
x = np.arange(0, Nx, 1)
y = np.arange(0, Ny, 1)
################################################################



################################################################
###                        SCATTERER                        ####
################################################################

a = a

# initilizing the polar coordinates
r, phi1 = [initialize_field(Nz, Ny, Nx) for _ in range(2)]

# center of the scatterer
cx = Nx//2 + 0.5
cy = Ny//2 + 0.5

# converting to polar coordinates
ModuleGeometry.carToPolar(r, phi, Nz, Ny, Nx, cy1, cx1)

# scatterer particle
ModuleGeometry.circle(r, phi, er_tot, er2, a1, Nz, Ny, Nx)

# coordinates where traction vector is being calculated
R = np.array([2, 3, 4, 5])*a

X_polar1 = R[0] * np.cos(theta*np.pi/180) + Nx//2 + 0.5
Y_polar1 = R[0] * np.sin(theta*np.pi/180) + Ny//2 + 0.5

X_polar2 = R[1] * np.cos(theta*np.pi/180) + Nx//2 + 0.5
Y_polar2 = R[1] * np.sin(theta*np.pi/180) + Ny//2 + 0.5

X_polar3 = R[2] * np.cos(theta*np.pi/180) + Nx//2 + 0.5
Y_polar3 = R[2] * np.sin(theta*np.pi/180) + Ny//2 + 0.5

X_polar4 = R[3] * np.cos(theta*np.pi/180) + Nx//2 + 0.5
Y_polar4 = R[3] * np.sin(theta*np.pi/180) + Ny//2 + 0.5

################################################################






###############################################################################################################
########                                           SHARED LIBRARY                                   ###########
###############################################################################################################

# loading the shared file (c library)
path = os.getcwd()
myclib = CDLL(os.path.join(path, "LBM_Sequential.so"))

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
myclib.collStream.argtypes = [P4D, P4D, P4D, P4D, P4D, P4D, P3D, P3D, P3D, P3D, P3D, P3D, P3D, P3D, c_int, c_int, c_int, c_int]
myclib.collStream.restype  = None

myclib.collStreamForcingNode.argtypes = [P4D, P4D, P4D, P4D, P4D, P4D, P3D, P3D, P3D, P3D, P3D, P3D, P3D, P3D, c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int]
myclib.collStreamForcingNode.restype  = None

###############################################################################################################



Ez_scat_interp1 = []
Ez_scat_interp2 = []
Ez_scat_interp3 = []
Ez_scat_interp4 = []


for t in range(int(Time)):

    t1 = time.time()



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
        planeWaveTM(Ez_inc, Hy_inc, t, omega, xmin, xmax, ymin, ymax)
        planeWaveTM(Ez_tot, Hy_tot, t, omega, xmin, xmax, ymin, ymax)

        # calculation of scattered fields
        Ez_scat = Ez_tot - Ez_inc
        Hx_scat = Hx_tot - Hx_inc
        Hy_scat = Hy_tot - Hy_inc


        # collision and streaming (the 2 steps of LBM) when field is forced
        myclib.collStreamForcingNode(ex_inc, ey_inc, ez_inc, hx_inc, hy_inc, hz_inc, Ex_inc, Ey_inc, Ez_inc, Hx_inc, Hy_inc, Hz_inc, er_inc, mur_inc, Nz, Ny, Nx, Q, xmin, xmax, ymin, ymax)
        myclib.collStreamForcingNode(ex_tot, ey_tot, ez_tot, hx_tot, hy_tot, hz_tot, Ex_tot, Ey_tot, Ez_tot, Hx_tot, Hy_tot, Hz_tot, er_tot, mur_tot, Nz, Ny, Nx, Q, xmin, xmax, ymin, ymax)

    else:

        # calculation of scattered fields
        Ez_scat = Ez_tot - Ez_inc
        Hx_scat = Hx_tot - Hx_inc
        Hy_scat = Hy_tot - Hy_inc

        
        # collision and streaming (the 2 steps of LBM) when field is not forced
        myclib.collStream(ex_inc, ey_inc, ez_inc, hx_inc, hy_inc, hz_inc, Ex_inc, Ey_inc, Ez_inc, Hx_inc, Hy_inc, Hz_inc, er_inc, mur_inc, Nz, Ny, Nx, Q)
        myclib.collStream(ex_tot, ey_tot, ez_tot, hx_tot, hy_tot, hz_tot, Ex_tot, Ey_tot, Ez_tot, Hx_tot, Hy_tot, Hz_tot, er_tot, mur_tot, Nz, Ny, Nx, Q)
        
    ###############################################################################################################



   
    

    if (t >= int(Time) - 5*int(np.round(period))):

        ###############################################################################################################
        ############################### Most time consuming part, 2/3rd of total time #################################
        
        # interpolating the field values
        Ez_scat_spline = RectBivariateSpline(x, y, Ez_scat[0])


        ###############################################################################################################
        ###############################################################################################################

        
        ###############################################################################################################
        ################################### This portion takes negligible time ########################################

        # interpolated values of fields at a circle of radius R
        Ez_scat_interp1.append(Ez_scat_spline.ev(Y_polar1, X_polar1))
        Ez_scat_interp2.append(Ez_scat_spline.ev(Y_polar2, X_polar2))
        Ez_scat_interp3.append(Ez_scat_spline.ev(Y_polar3, X_polar3))
        Ez_scat_interp4.append(Ez_scat_spline.ev(Y_polar4, X_polar4))

        ###############################################################################################################
        ###############################################################################################################


    if (t == 10):
        t2 = time.time()
        print('Approximate time of computation in minutes :', int((t2 - t1)*Time/60))
    


    
##############################################################

EzScat1 = open(directory+"/Ez_scat1.txt", "w")
np.savetxt(EzScat1, Ez_scat_interp1, fmt='%.4e')
EzScat1.close()

EzScat2 = open(directory+"/Ez_scat2.txt", "w")
np.savetxt(EzScat2, Ez_scat_interp2, fmt='%.4e')
EzScat2.close()

EzScat3 = open(directory+"/Ez_scat3.txt", "w")
np.savetxt(EzScat3, Ez_scat_interp3, fmt='%.4e')
EzScat3.close()

EzScat4 = open(directory+"/Ez_scat4.txt", "w")
np.savetxt(EzScat4, Ez_scat_interp4, fmt='%.4e')
EzScat4.close()


    

