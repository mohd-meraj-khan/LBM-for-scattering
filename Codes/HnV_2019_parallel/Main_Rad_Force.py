from ctypes import CDLL, POINTER
from ctypes import c_size_t, c_double, c_float, c_int
import numpy as np
from scipy.interpolate import RectBivariateSpline
import sys
import os
import time

from Module_Traction import tractionFx, tractionFy, tractionFz
from Module_EM_Wave import planeWaveTM
import Module_Geometry

from parameters import *
from Shared_Lib import *


t0 = time.time()



directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)



######################
# LBM properties (DO NOT CHANGE)
Q = 12               # number of velocities at a grid
velocity = 1.0/2    # velocity of EM wave in vacuum
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

Px, Py, Pz, Mx, My, Mz = [initialize_field(Ny, Nx) for _ in range(6)]

Pxb, Pyb, Pzb, Mxb, Myb, Mzb = [initialize_field(Ny, Nx) for _ in range(6)]


# initializing the distribution functions of electric and magnetic fields
def initilize_dis_func(Ny=10, Nx=10, Q=12):
    return np.zeros((Ny, Nx, Q), dtype=np.float32, order='C')

f, fb = [initilize_dis_func(Ny, Nx, Q) for _ in range(2)]



# initilizing the domain properties
def initialize_material_properties(er1=1, mur1=1, Ny=10, Nx=10):
    yield np.ones((Ny, Nx), dtype=np.float32, order='C') * er1
    yield np.ones((Ny, Nx), dtype=np.float32, order='C') * mur1
    return

er, mur = initialize_material_properties(er1, mur1, Ny, Nx)
#############################################################################################################




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

# coordinates where traction vector is being calculated
R = 1.0*a
X_polar = R * np.cos(theta*np.pi/180) + cx
Y_polar = R * np.sin(theta*np.pi/180) + cy

################################################################





t1 = time.time()



Fx1 = []
Fx_avg1 = []
Traction1 = []



for t in range(int(Time)):
    
    

    #################################################################################################################
    ########                                         LBM CALCULATION                                          #######
    #################################################################################################################

    # computation of macroscopic fields from distribution function
    myclib.macroField(f, Ex, Ey, Ez, Hx, Hy, Hz, Px, Py, Pz, Mx, My, Mz, er, mur, Ny, Nx, Q, N)


    
    if (t >= 0):
        
        # source wave
        planeWaveTM(Ez, Hy, t, omega, xloc, ymin, ymax)

        # collision and streaming (the 2 steps of LBM) when field is forced
        myclib.collForcingNode(f, fb, Ex, Ey, Ez, Hx, Hy, Hz, Px, Py, Pz, Mx, My, Mz, Pxb, Pyb, Pzb, Mxb, Myb, Mzb, er, mur, Ny, Nx, Q, xloc, ymin, ymax, N)
        myclib.streaming(f, fb, Px, Py, Pz, Mx, My, Mz, Pxb, Pyb, Pzb, Mxb, Myb, Mzb, Ny, Nx, Q, N)
        

    else:
        
        # collision and streaming (the 2 steps of LBM) when field is not forced
        myclib.collNotForcingNode(f, fb, Ex, Ey, Ez, Hx, Hy, Hz, Px, Py, Pz, Mx, My, Mz, Pxb, Pyb, Pzb, Mxb, Myb, Mzb, er, mur, Ny, Nx, Q, N)
        myclib.streaming(f, fb, Px, Py, Pz, Mx, My, Mz, Pxb, Pyb, Pzb, Mxb, Myb, Mzb, Ny, Nx, Q, N)
        
    ###############################################################################################################

   
    
    if (t >= int(Time) - int(np.round(period))):

        ###############################################################################################################
        
        # interpolating the field values
        Ex_spline = RectBivariateSpline(x, y, Ex)
        Ey_spline = RectBivariateSpline(x, y, Ey)
        Ez_spline = RectBivariateSpline(x, y, Ez)

        Hx_spline = RectBivariateSpline(x, y, Hx)
        Hy_spline = RectBivariateSpline(x, y, Hy)
        Hz_spline = RectBivariateSpline(x, y, Hz)
        ###############################################################################################################
        


        ###############################################################################################################
        ###############                         FORCE CALCULATION AT SCATTERER 1                        ###############
        ###############################################################################################################

        # interpolated values of fields at a circle of radius R1
        Ex_interp = Ex_spline.ev(Y_polar, X_polar)
        Ey_interp = Ey_spline.ev(Y_polar, X_polar)
        Ez_interp = Ez_spline.ev(Y_polar, X_polar)

        Hx_interp = Hx_spline.ev(Y_polar, X_polar)
        Hy_interp = Hy_spline.ev(Y_polar, X_polar)
        Hz_interp = Hz_spline.ev(Y_polar, X_polar)
        
        # integrating the traction force over the perimeter of the circle
        TdotN1 = tractionFx(Ex_interp, Ey_interp, Ez_interp, Hx_interp, Hy_interp, Hz_interp, er1, mur1, nx, ny, nz)
        Traction1.append(TdotN1)
        Fx1.append(np.sum(TdotN1 * R * dtheta * (np.pi / 180) / (a / ratio)))
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
