import pyvista as pv
import numpy as np
import os
import time

from Module_Traction_3D import *
from Module_EM_Wave_3D import *
from Module_Parameters_3D import *
from Module_Shared_Lib_3D import *

t0 = time.time()

pictures = 'gallery'
if not os.path.exists(pictures):
    os.makedirs(pictures)

print("Number of time steps :", int(Time))

###############################################################################################################
####      DEFINING AND INITILIZING VARIABLES FOR MACROSCOPIC FIELDS AND DISTRIBUTION FUNCTIONS             ####
###############################################################################################################

'''initializing the electric and magnetic fields'''
def initialize_field(Nz=10, Ny=10, Nx=10):
    return np.zeros((Nz, Ny, Nx), dtype=np.float32, order='C')

Ex, Ey, Ez, Hx, Hy, Hz = [initialize_field(Nz, Ny, Nx) for _ in range(6)]

'''initializing the distribution functions of electric and magnetic fields'''
def initilize_dis_func(Nz=10, Ny=10, Nx=10, Q=7):
    return np.zeros((Nz, Ny, Nx, Q), dtype=np.float32, order='C')

ex, ey, ez, hx, hy, hz = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]
exb, eyb, ezb, hxb, hyb, hzb = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]

###############################################################################################################

fps = 30
Video_Name = "Ez_tot.mp4"

x = np.arange(Nx, dtype=np.float32)
y = np.arange(Ny, dtype=np.float32)
z = np.arange(Nz, dtype=np.float32)
Z, Y, X = np.meshgrid(z, y, x, indexing="ij")

points = np.vstack((X.ravel(), Y.ravel(), Z.ravel())).T

# Create PolyData and Plotter once
poly_data = pv.PolyData(points)
poly_data.point_data["Ez"] = Ez.ravel()

plotter = pv.Plotter()
plotter.show_axes()
plotter.open_movie(Video_Name, framerate=fps)

# Add the mesh to the plotter
plotter.add_mesh(poly_data, scalars="Ez", cmap="coolwarm", point_size=10)

# Update field function for time-stepping
def update_field(t):
    """ Update the field data and visualization at each time step """
    myclib.initializeField(Ex, Ey, Ez, Hx, Hy, Hz, Nz, Ny, Nx, N)

    myclib.macroField(ex, er, Ex, Nz, Ny, Nx, Q, N)
    myclib.macroField(ey, er, Ey, Nz, Ny, Nx, Q, N)
    myclib.macroField(ez, er, Ez, Nz, Ny, Nx, Q, N)

    myclib.macroField(hx, mur, Hx, Nz, Ny, Nx, Q, N)
    myclib.macroField(hy, mur, Hy, Nz, Ny, Nx, Q, N)
    myclib.macroField(hz, mur, Hz, Nz, Ny, Nx, Q, N)

    if t >= 0:
        planeWaveTM(Ez, Hy, t, omega, xloc, ymin, ymax, zmin, zmax)
        myclib.collForcingNode(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Nz, Ny, Nx, Q, xloc, ymin, ymax, zmin, zmax, N)
        myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Nz, Ny, Nx, Q, N)
    else:
        myclib.collision(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Nz, Ny, Nx, Q, N)
        myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Nz, Ny, Nx, Q, N)

    # Update the scalar data
    poly_data.point_data["Ez"] = Ez.ravel()

# Run the simulation and animation
t1 = time.time()

for t in range(int(Time)):
    update_field(t)
    plotter.write_frame()  # Write the current frame

    # Time estimation
    t2 = time.time()
    if t > 0 and t % 100 == 0:
        remaining_time = (t2 - t1) * (int(Time) - t) / (t * 60)
        print(f"Approximate time left: {remaining_time:.2f} minutes", end="\r")

# Close the plotter to finalize the video
plotter.close()

t3 = time.time()
print(f"\nTotal time taken: {(t3 - t0) / 60:.2f} minutes")
