import pyvista as pv

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


##pv.OFF_SCREEN = True



fps = 30

Video_Name = "Ez_tot.mp4"


x = np.arange(Nx, dtype=np.float32)
y = np.arange(Ny, dtype=np.float32)
z = np.arange(Nz, dtype=np.float32)
Z, Y, X = np.meshgrid(z, y, x, indexing="ij")

# PyVista structured grid
grid = pv.StructuredGrid(X, Y, Z)
grid.point_data['Ez'] = Ez.ravel()  # Attach Ez field data to the grid

# Set up the plotter and open a movie with desired FPS
plotter = pv.Plotter()
actor = plotter.add_mesh(grid, scalars='Ez', colormap='seismic', opacity=1.0, clim=[-1, 1])
plotter.show_axes()
plotter.open_movie(Video_Name)  # Set FPS


# Show axes
plotter.show_axes()

# Open movie file
plotter.open_movie(Video_Name, framerate=fps)


    
        
def update_field(t):   

    #################################################################################################################
    ########                                         LBM CALCULATION                                          #######
    #################################################################################################################

    '''initialization of macroscopic fields'''
    myclib.initializeField(Ex, Ey, Ez, Hx, Hy, Hz, Nz, Ny, Nx, N)

    '''computation of macroscopic fields from distribution function'''
    myclib.macroField(ex, er, Ex, Nz, Ny, Nx, Q, N)
    myclib.macroField(ey, er, Ey, Nz, Ny, Nx, Q, N)
    myclib.macroField(ez, er, Ez, Nz, Ny, Nx, Q, N)

    myclib.macroField(hx, mur, Hx, Nz, Ny, Nx, Q, N)
    myclib.macroField(hy, mur, Hy, Nz, Ny, Nx, Q, N)
    myclib.macroField(hz, mur, Hz, Nz, Ny, Nx, Q, N)

    
    if (t >= 0):
        
        '''source wave'''
        planeWaveTM(Ez, Hy, t, omega, xloc, ymin, ymax, zmin, zmax)

        '''collision and streaming (the 2 steps of LBM) when field is forced'''
        myclib.collForcingNode(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Nz, Ny, Nx, Q, xloc, ymin, ymax, zmin, zmax, N)
        myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Nz, Ny, Nx, Q, N)
    else:    
        '''collision and streaming (the 2 steps of LBM) when field is not forced'''
        myclib.collision(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Nz, Ny, Nx, Q, N)
        myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Nz, Ny, Nx, Q, N)
            
    ###############################################################################################################

    grid.point_data['Ez'] = Ez.ravel()
    actor.mapper.update()
    # Extract three perpendicular slices using slice_orthogonal
    slices = grid.slice_orthogonal(x=50, y=50, z=50)

    # Clear previous plot and add the orthogonal slices
    plotter.clear()
    plotter.add_mesh(slices, scalars='Ez', colormap='seismic', clim=[-1, 1])
    plotter.camera.azimuth += 0.1
    plotter.camera.elevation = np.sin(t * 0.01) * 50
##    plotter.camera.elevation += 0.5

t1 = time.time()


for t in range(int(Time)):

    


    update_field(t)
    plotter.write_frame()  # Write the current frame to the movie




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
