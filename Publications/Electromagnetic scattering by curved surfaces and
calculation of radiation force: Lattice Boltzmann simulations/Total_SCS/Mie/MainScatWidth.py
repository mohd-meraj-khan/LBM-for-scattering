from ctypes import CDLL, POINTER
from ctypes import c_size_t, c_double, c_float, c_int
import numpy as np
from scipy.interpolate import RectBivariateSpline
import sys
import os
import time


from ModuleSourceOfEmWave import planeWaveTM
import ModuleGeometry


t0 = time.time()



directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)




# Accessing command line arguments
parameters = sys.argv

ratio = float(parameters[1])

######################


a = 100
n = 4
Nx, Ny, Nz = n*a, n*a, 1  # size of the computational domain

er1, mur1, er2 = 1, 1, 4   # material properties i.e. permittivity and permeabilty


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
noOfReflections = 200

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
################### for interpolation ##########################

dtheta = 1
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
r, phi = [initialize_field(Nz, Ny, Nx) for _ in range(2)]

# center of the scatterer
cx = Nx//2 + 0.5
cy = Ny//2 + 0.5

# converting to polar coordinates
ModuleGeometry.carToPolar(r, phi, Nz, Ny, Nx, cy, cx)

# scatterer particle
scatterer = ModuleGeometry.circle(r, a, Nz, Ny, Nx)

er_tot[scatterer] = er2

# coordinates where traction vector is being calculated
R = 1.0*a

X_polar = R * np.cos(theta*np.pi/180) + cx
Y_polar = R * np.sin(theta*np.pi/180) + cy



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


Ez_scat_interp = []
Hx_scat_interp = []
Hy_scat_interp = []



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

        
        # collision and streaming (the 2 steps of LBM) when field is not forced
        myclib.collStream(ex_inc, ey_inc, ez_inc, hx_inc, hy_inc, hz_inc, exb_inc, eyb_inc, ezb_inc, hxb_inc, hyb_inc, hzb_inc, Ex_inc, Ey_inc, Ez_inc, Hx_inc, Hy_inc, Hz_inc, er_inc, mur_inc, Nz, Ny, Nx, Q)
        myclib.collStream(ex_tot, ey_tot, ez_tot, hx_tot, hy_tot, hz_tot, exb_tot, eyb_tot, ezb_tot, hxb_tot, hyb_tot, hzb_tot, Ex_tot, Ey_tot, Ez_tot, Hx_tot, Hy_tot, Hz_tot, er_tot, mur_tot, Nz, Ny, Nx, Q)
        myclib.streaming(ex_inc, ey_inc, ez_inc, hx_inc, hy_inc, hz_inc, exb_inc, eyb_inc, ezb_inc, hxb_inc, hyb_inc, hzb_inc, Nz, Ny, Nx, Q)
        myclib.streaming(ex_tot, ey_tot, ez_tot, hx_tot, hy_tot, hz_tot, exb_tot, eyb_tot, ezb_tot, hxb_tot, hyb_tot, hzb_tot, Nz, Ny, Nx, Q)
        
    ###############################################################################################################

   
    

    if (t >= int(Time) - 5*int(np.round(period))):

        ###############################################################################################################
        
        # interpolating the field values
        Ez_scat_spline = RectBivariateSpline(x, y, Ez_scat[0])
        Hx_scat_spline = RectBivariateSpline(x, y, Hx_scat[0])
        Hy_scat_spline = RectBivariateSpline(x, y, Hy_scat[0])



        # interpolated values of fields at a circle of radius R
        Ez_scat_interp.append(Ez_scat_spline.ev(Y_polar, X_polar))
        Hx_scat_interp.append(Hx_scat_spline.ev(Y_polar, X_polar))
        Hy_scat_interp.append(Hy_scat_spline.ev(Y_polar, X_polar))
        

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
    


    
##############################################################

EzScat = open(directory+"/Ez_scat_{}.txt".format(ratio), "w")
np.savetxt(EzScat, Ez_scat_interp, fmt='%.4e')
EzScat.close()


HxScat = open(directory+"/Hx_scat_{}.txt".format(ratio), "w")
np.savetxt(HxScat, Hx_scat_interp, fmt='%.4e')
HxScat.close()


HyScat = open(directory+"/Hy_scat_{}.txt".format(ratio), "w")
np.savetxt(HyScat, Hy_scat_interp, fmt='%.4e')
HyScat.close()
    

