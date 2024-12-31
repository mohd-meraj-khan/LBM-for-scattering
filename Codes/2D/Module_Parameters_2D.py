import sys
import time
import numpy as np
from Module_Geometry_2D import *



'''number of parallel threads'''
N = 12


'''Accessing command line arguments'''
parameters = sys.argv

theta = 45
######################

a, ratio = 25, 0.5   # ratio = a / wavelength
n = 6
Nx, Ny = n*a, n*a  # size of the computational domain

er1, mur1, er2, er3 = 1, 1, 4, 10000   # material properties i.e. permittivity and permeabilty


################################################################
###                        SCATTERER                        ####
################################################################

'''initilizing the polar coordinates'''
def initialize_field(Ny=10, Nx=10):
    return np.zeros((Ny, Nx), dtype=np.float32, order='C')
r, phi = [initialize_field(Ny, Nx) for _ in range(2)]


'''initilizing the domain properties'''
def initialize_material_properties(er1=1, mur1=1, Ny=10, Nx=10):
    yield np.ones((Ny, Nx), dtype=np.float32, order='C') * er1
    yield np.ones((Ny, Nx), dtype=np.float32, order='C') * mur1
    return

er, mur = initialize_material_properties(er1, mur1, Ny, Nx)

'''center of the scatterer'''
cx = Nx//2 + 0.5
cy = Ny//2 + 0.5

'''converting from cartesian to polar coordinates'''
carToPolar(r, phi, Ny, Nx, cy, cx)

'''scatterer particle'''
##first_half, second_half  = Module_Geometry.JanusHalf(r, a, theta, Ny, Nx, cy, cx)
##
##er[first_half]  = er2
##er[second_half] = er3


scatterer = circle(r, a, Ny, Nx)

er[scatterer] = er2

##er = (er2 - er1)/2 * (np.tanh(a - r) + 1) + er1
##
##
##print(er[Nx//2, :])
##sys.exit()
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
#############################



'''half-width of the square bounding box'''
w = int(np.round(1.5*a))

'''box surrounding the scatterer'''
Top    = int(cy + w), np.arange(int(cx - w), int(cx + w))
Right  = np.arange(int(cy - w), int(cy + w)), int(cx + w)
Bottom = int(cy - w), np.arange(int(cx - w), int(cx + w))
Left   = np.arange(int(cy - w), int(cy + w)), int(cx - w)




'''unit normal vectors at the perimeter of the bounding box'''
nxTop, nxRight, nxBottom, nxLeft = 0, 1, 0, -1
nyTop, nyRight, nyBottom, nyLeft = 1, 0, -1, 0
nzTop, nzRight, nzBottom, nzLeft = 0, 0,  0, 0







noOfPeriods = 0
noOfReflections = 10


'''number of time steps the code will run'''
Time = int(3 * (Nx * np.sqrt(er1) + noOfReflections * 2 * a * np.sqrt(er2)) + noOfPeriods * period)

