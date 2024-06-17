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

ratio = float(parameters[1])

######################

#ratio = 0.95   # ratio = a / wavelength

a  = 25
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

# coordinates where traction vector is being calculated
R1 = 1.0*a1
X_polar1 = R1 * np.cos(theta*np.pi/180) + cx1
Y_polar1 = R1 * np.sin(theta*np.pi/180) + cy1

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






Fx1 = []
Fx_avg1 = []
Traction1 = []



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

    # plane wave enforcement at the left boundary of the domain
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

    t2 = time.time()



    
    if (t >= int(Time) - int(np.round(period))):

        ###############################################################################################################
        
        # interpolating the field values
        Ex_spline = RectBivariateSpline(x, y, Ex[0])
        Ey_spline = RectBivariateSpline(x, y, Ey[0])
        Ez_spline = RectBivariateSpline(x, y, Ez[0])

        Hx_spline = RectBivariateSpline(x, y, Hx[0])
        Hy_spline = RectBivariateSpline(x, y, Hy[0])
        Hz_spline = RectBivariateSpline(x, y, Hz[0])
        ###############################################################################################################
        


        ###############################################################################################################
        ###############                         FORCE CALCULATION AT SCATTERER 1                        ###############
        ###############################################################################################################

        # interpolated values of fields at a circle of radius R1
        Ex_interp = Ex_spline.ev(Y_polar1, X_polar1)
        Ey_interp = Ey_spline.ev(Y_polar1, X_polar1)
        Ez_interp = Ez_spline.ev(Y_polar1, X_polar1)

        Hx_interp = Hx_spline.ev(Y_polar1, X_polar1)
        Hy_interp = Hy_spline.ev(Y_polar1, X_polar1)
        Hz_interp = Hz_spline.ev(Y_polar1, X_polar1)
        
        # integrating the traction force over the perimeter of the circle
        TdotN1 = tractionFx(Ex_interp, Ey_interp, Ez_interp, Hx_interp, Hy_interp, Hz_interp, er1, mur1, nx, ny, nz)
        Traction1.append(TdotN1)
        Fx1.append(np.sum(TdotN1 * R1 * dtheta * (np.pi / 180) / (a1 / ratio)))
        ###############################################################################################################
    



Fx_avg1.append(np.average(Fx1))
print(Fx_avg1[0])



tracVec1 = open(directory+"/tracVec1_ratio_{}.txt".format(ratio), "w")
np.savetxt(tracVec1, Traction1, fmt='%.4e')
tracVec1.close()

fxSum1 = open(directory+"/fxSum1_ratio_{}.txt".format(ratio), "w")
np.savetxt(fxSum1, Fx1, fmt='%.4e')
fxSum1.close()

fxAvg1 = open(directory+"/FxAvg1.txt", "a")
np.savetxt(fxAvg1, Fx_avg1, fmt='%.4e')
fxAvg1.close()
