import sys
import time
import numpy as np
from Module_Geometry_3D import *

import time

'''number of parallel threads'''
N = 10


'''Accessing command line arguments'''
parameters = sys.argv


phi0 = 0

######################

ratio = 5
er2 = 2


er1, mur1 = 1, 1   # material properties i.e. permittivity and permeabilty

mur2 = 1

V1 = 1 / (3 * np.sqrt(er1*mur1))
V2 = 1 / (3 * np.sqrt(er2*mur2))


A = 30

if (ratio <= 1 * V2 / V1):
    a = A
else:
    a = int(np.round(A * ratio * V1 / V2))




if (ratio < 1):
    n = 10
else:
    n = 3


Nx, Ny, Nz = n*a, n*a, n*a  # size of the computational domain



################################################################
###                        SCATTERER                        ####
################################################################

'''initilizing the polar coordinates'''
def initialize_field(Nz=10, Ny=10, Nx=10):
    return np.zeros((Nz, Ny, Nx), dtype=np.float32, order='C')
rad, the, phi = [initialize_field(Nz, Ny, Nx) for _ in range(3)]


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
rad, the, phi = carToSpherical(Nz, Ny, Nx, cz, cy, cx)

'''scatterer particle'''
scatterer = sphere(rad, a)

er[scatterer] = er2



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
zloc = 0
ymin = 0
ymax = Ny
xmin = 0
xmax = Nx
#############################



'''half-width of the square bounding box'''
w = int(np.round(1.25*a))

'''box surrounding the scatterer'''
Top = (int(cz + w), slice(int(cy - w), int(cy + w)), slice(int(cx - w), int(cx + w)))
Bot = (int(cz - w), slice(int(cy - w), int(cy + w)), slice(int(cx - w), int(cx + w)))
Ryt = (slice(int(cz - w), int(cz + w)), int(cy + w), slice(int(cx - w), int(cx + w)))
Lef = (slice(int(cz - w), int(cz + w)), int(cy - w), slice(int(cx - w), int(cx + w)))
Frt = (slice(int(cz - w), int(cz + w)), slice(int(cy - w), int(cy + w)), int(cx + w))
Bak = (slice(int(cz - w), int(cz + w)), slice(int(cy - w), int(cy + w)), int(cx - w))


'''unit normal vectors at the surface of the bounding box'''
nxTop, nxBot, nxRyt, nxLef, nxFrt, nxBak = 0, 0, 0, 0, 1, -1
nyTop, nyBot, nyRyt, nyLef, nyFrt, nyBak = 0, 0, 1, -1, 0, 0
nzTop, nzBot, nzRyt, nzLef, nzFrt, nzBak = 1, -1, 0, 0, 0, 0



'''x, y and z coordinates at the surface of the bounding box to compute torque'''
x_id_top = np.arange(int(cx-w), int(cx+w))
x_id_bot = np.arange(int(cx-w), int(cx+w))
y_id_top = np.arange(int(cy-w), int(cy+w))
y_id_bot = np.arange(int(cy-w), int(cy+w))

x_id_ryt = np.arange(int(cx-w), int(cx+w))
x_id_lef = np.arange(int(cx-w), int(cx+w))
z_id_ryt = np.arange(int(cz-w), int(cz+w))
z_id_lef = np.arange(int(cz-w), int(cz+w))

y_id_frt = np.arange(int(cy-w), int(cy+w))
y_id_bak = np.arange(int(cy-w), int(cy+w))
z_id_frt = np.arange(int(cz-w), int(cz+w))
z_id_bak = np.arange(int(cz-w), int(cz+w))

X_top, Y_top = np.meshgrid(x_id_top, y_id_top, indexing='xy')
X_bot, Y_bot = np.meshgrid(x_id_bot, y_id_bot, indexing='xy')

Z_top = int(cz + w)
Z_bot = int(cz - w)

X_ryt, Z_ryt = np.meshgrid(x_id_ryt, z_id_ryt, indexing='xy')
X_lef, Z_lef = np.meshgrid(x_id_lef, z_id_lef, indexing='xy')

Y_ryt = int(cy + w)
Y_lef = int(cy - w)

Y_frt, Z_frt = np.meshgrid(y_id_frt, z_id_frt, indexing='xy')
Y_bak, Z_bak = np.meshgrid(y_id_bak, z_id_bak, indexing='xy')

X_frt = int(cx + w)
X_bak = int(cx - w)









if (ratio < 1):
    noOfPeriods = 5
else:
    noOfPeriods = 20



noOfReflections = 0


'''number of time steps the code will run'''
Time = int(3 * (Nx * np.sqrt(er1) + noOfReflections * 2 * a * np.sqrt(er2)) + noOfPeriods * period)

