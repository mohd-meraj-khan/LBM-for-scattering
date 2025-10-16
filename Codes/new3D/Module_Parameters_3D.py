import sys
import time
import numpy as np
from Module_Geometry_3D import *



'''number of parallel threads'''
N = 4


'''Accessing command line arguments'''
parameters = sys.argv

theta = 45
######################

a, ratio = 25, 0.5   # ratio = a / wavelength
n = 4
Nx, Ny, Nz = n*a, n*a, n*a  # size of the computational domain

er1, mur1, er2, er3 = 1, 1, 36, 10000   # material properties i.e. permittivity and permeabilty


################################################################
###                        SCATTERER                        ####
################################################################

'''initilizing the polar coordinates'''
def initialize_field(Nz=10, Ny=10, Nx=10):
    return np.zeros((Nz, Ny, Nx), dtype=np.float32, order='C')
r, theta, phi = [initialize_field(Nz, Ny, Nx) for _ in range(3)]


'''initilizing the domain properties'''
def initialize_material_properties(er1=1, mur1=1, Nz=10, Ny=10, Nx=10):
    yield np.ones((Nz, Ny, Nx), dtype=np.float32, order='C') * er1
    yield np.ones((Nz, Ny, Nx), dtype=np.float32, order='C') * mur1
    return

er, mur = initialize_material_properties(er1, mur1, Nz, Ny, Nx)

'''center of the scatterer'''
cx = Nx//2 + 0.5
cy = Ny//2 + 0.5
cz = Nz//2 + 0.5

'''converting from cartesian to spherical coordinates'''
r, theta, phi = carToSpherical(Nz, Ny, Nx, cz, cy, cx)

'''scatterer particle'''
sphere(r, er, er2, a)



################################################################


######################
# LBM properties (DO NOT CHANGE)
Q = 7               # number of velocities at a grid
velocity = 1.0/3    # velocity of EM wave in vacuum
######################


######### DON'T CHANGE ########
wavelength = a / ratio
omega = 2 * np.pi * velocity / wavelength
period = int(np.round(wavelength / velocity))

'''wavenumber of the incident wave'''
k = 2*np.pi / wavelength
###############################


#############################
# boundary of EM wave source
xloc = 0
ymin = 0
ymax = Ny
zmin = 0
zmax = Nz
#############################



'''half-width of the square bounding box'''
w = int(np.round(1.5*a))

'''box surrounding the scatterer'''
Top    = np.arange(int(cz - w), int(cz + w)), int(cy + w), np.arange(int(cx - w), int(cx + w))
Bottom = np.arange(int(cz - w), int(cz + w)), int(cy - w), np.arange(int(cx - w), int(cx + w))
Right  = np.arange(int(cz - w), int(cz + w)), np.arange(int(cy - w), int(cy + w)), int(cx + w)
Left   = np.arange(int(cz - w), int(cz + w)), np.arange(int(cy - w), int(cy + w)), int(cx - w)
Front  = int(cz + w), np.arange(int(cy - w), int(cy + w)), np.arange(int(cx - w), int(cx + w))
Back   = int(cz - w), np.arange(int(cy - w), int(cy + w)), np.arange(int(cx - w), int(cx + w))

'''unit normal vectors at the perimeter of the bounding box'''
nxTop, nxBottom, nxRight, nxLeft, nxFront, nxBack = 0, 0, 1, -1, 0, 0
nyTop, nyBottom, nyRight, nyLeft, nyFront, nyBack = 1, -1, 0, 0, 0, 0
nzTop, nzBottom, nzRight, nzLeft, nzFront, nzBack = 0, 0, 0, 0, 1, -1







noOfPeriods = 10
noOfReflections = 0


'''number of time steps the code will run'''
Time = int(3 * (Nx * np.sqrt(er1) + noOfReflections * 2 * a * np.sqrt(er2)) + noOfPeriods * period)

