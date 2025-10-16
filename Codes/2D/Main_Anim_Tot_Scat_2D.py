import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec

from Module_Traction_2D import *
from Module_EM_Wave_2D import *
from Module_Parameters_2D import *
from Module_Shared_Lib_2D import *


t0 = time.time()


pictures = 'gallery'
if not os.path.exists(pictures):
    os.makedirs(pictures)




print('\n')
print(f'ratio : {ratio}')
print(f'Number of parallel threads :{N}')
print(f'Number of time steps :{Time}')

print(f"Size of the computational domain: {Ny} * {Nx}\n\n")

print(f"Radius of the cylinder: {a}.\n")
print(f"Wavelength of the incident wave: {wavelength:.2f}.\n\n")





###############################################################################################################
####      DEFINING AND INITILIZING VARIABLES FOR MACROSCOPIC FIELDS AND DISTRIBUTION FUNCTIONS             ####
###############################################################################################################

'''initializing the electric and magnetic fields'''
def initialize_field(Ny=10, Nx=10):
    return np.zeros((Ny, Nx), dtype=np.float32, order='C')

ExI, EyI, EzI, HxI, HyI, HzI = [initialize_field(Ny, Nx) for _ in range(6)]
Ex, Ey, Ez, Hx, Hy, Hz = [initialize_field(Ny, Nx) for _ in range(6)]
Ex_scat, Ey_scat, Ez_scat, Hx_scat, Hy_scat, Hz_scat = [initialize_field(Ny, Nx) for _ in range(6)]


'''initializing the distribution functions of electric and magnetic fields'''
def initilize_dis_func(Ny=10, Nx=10, Q=7):
    return np.zeros((Ny, Nx, Q), dtype=np.float32, order='C')

exI, eyI, ezI, hxI, hyI, hzI = [initilize_dis_func(Ny, Nx, Q) for _ in range(6)]
exbI, eybI, ezbI, hxbI, hybI, hzbI = [initilize_dis_func(Ny, Nx, Q) for _ in range(6)]
ex, ey, ez, hx, hy, hz = [initilize_dis_func(Ny, Nx, Q) for _ in range(6)]
exb, eyb, ezb, hxb, hyb, hzb = [initilize_dis_func(Ny, Nx, Q) for _ in range(6)]

'''initializing the material properties for incident fields'''
erI, murI = initialize_material_properties(er1, mur1, Ny, Nx)
###############################################################################################################







t1 = time.time()


fig = plt.figure(figsize = (3.35*3, 3.35), dpi=600)
gs=GridSpec(1,2)
plt.ion()

fps = 10

Video_Name = "EzTotScat.mp4"

FFMpegWriter = animation.writers['ffmpeg']

metadata = dict(title='Ez', artist='', comment='Movie support!')
writer = FFMpegWriter(fps=fps, metadata=metadata)

with writer.saving(fig, Video_Name,300):
    for t in range(int(Time)):
        
        

        #################################################################################################################
        ########                                         LBM CALCULATION                                          #######
        #################################################################################################################

        '''initialization of macroscopic fields'''
        myclib.initializeField(ExI, EyI, EzI, HxI, HyI, HzI, Ny, Nx, N)
        myclib.initializeField(Ex, Ey, Ez, Hx, Hy, Hz, Ny, Nx, N)

        '''computation of macroscopic fields from distribution function'''
        myclib.macroField(exI, erI, ExI, Ny, Nx, Q, N)
        myclib.macroField(eyI, erI, EyI, Ny, Nx, Q, N)
        myclib.macroField(ezI, erI, EzI, Ny, Nx, Q, N)

        myclib.macroField(hxI, murI, HxI, Ny, Nx, Q, N)
        myclib.macroField(hyI, murI, HyI, Ny, Nx, Q, N)
        myclib.macroField(hzI, murI, HzI, Ny, Nx, Q, N)


        myclib.macroField(ex, er, Ex, Ny, Nx, Q, N)
        myclib.macroField(ey, er, Ey, Ny, Nx, Q, N)
        myclib.macroField(ez, er, Ez, Ny, Nx, Q, N)

        myclib.macroField(hx, mur, Hx, Ny, Nx, Q, N)
        myclib.macroField(hy, mur, Hy, Ny, Nx, Q, N)
        myclib.macroField(hz, mur, Hz, Ny, Nx, Q, N)


        if (t <= period):
            
            '''source wave'''
            planeWaveTM(EzI, HyI, t, omega, xloc, ymin, ymax)
            planeWaveTM(Ez, Hy, t, omega, xloc, ymin, ymax)


            '''calculation of scattered fields'''
            Ex_scat = Ex - ExI
            Ey_scat = Ey - EyI
            Ez_scat = Ez - EzI
            Hx_scat = Hx - HxI
            Hy_scat = Hy - HyI
            Hz_scat = Hz - HzI

            '''electric field is zero inside PEC'''
##            Ex[scatterer] = 0
##            Ey[scatterer] = 0
##            Ez[scatterer] = 0

            '''collision and streaming (the 2 steps of LBM) when field is forced'''
            myclib.collForcingNode(exI, eyI, ezI, hxI, hyI, hzI, exbI, eybI, ezbI, hxbI, hybI, hzbI, ExI, EyI, EzI, HxI, HyI, HzI, erI, murI, Ny, Nx, Q, xloc, ymin, ymax, N)
            myclib.collForcingNode(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Ny, Nx, Q, xloc, ymin, ymax, N)
            myclib.streaming(exI, eyI, ezI, hxI, hyI, hzI, exbI, eybI, ezbI, hxbI, hybI, hzbI, Ny, Nx, Q, N)
            myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ny, Nx, Q, N)



        else:
            '''calculation of scattered fields'''
            Ex_scat = Ex - ExI
            Ey_scat = Ey - EyI
            Ez_scat = Ez - EzI
            Hx_scat = Hx - HxI
            Hy_scat = Hy - HyI
            Hz_scat = Hz - HzI


            '''electric field is zero inside PEC'''
##            Ex[scatterer] = 0
##            Ey[scatterer] = 0
##            Ez[scatterer] = 0

            
            '''collision and streaming (the 2 steps of LBM) when field is not forced'''
            myclib.collNotForcingNode(exI, eyI, ezI, hxI, hyI, hzI, exbI, eybI, ezbI, hxbI, hybI, hzbI, ExI, EyI, EzI, HxI, HyI, HzI, erI, murI, Ny, Nx, Q, N)
            myclib.collNotForcingNode(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Ny, Nx, Q, N)
            myclib.streaming(exI, eyI, ezI, hxI, hyI, hzI, exbI, eybI, ezbI, hxbI, hybI, hzbI, Ny, Nx, Q, N)
            myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ny, Nx, Q, N)
        ###############################################################################################################




        

    ###############################################################################################################
    ##########                                           ANIMATION                                       ##########
    ###############################################################################################################

        if (t%10 == 0):
                
            fig.clear()

            ########### Plot parameters ################
            plt.rc('font', family = 'serif', size = 10)
            plt.rc('xtick', labelsize = 10)
            plt.rc('ytick', labelsize = 10)
            plt.rc('lines', markersize = 2, lw = 0.75)
            plt.rc('text', usetex = True)
            plt.rcParams['image.cmap']='seismic'
            ###########################################
            
            
                        
            ax1 = fig.add_subplot(gs[0,0])
            plt.title(r'$E_z^{tot}$')
            im1 = plt.imshow(Ez, vmin = -1, vmax = 1, cmap='seismic', origin='lower')
            im = plt.imshow(er, extent=(0, Nx, 0, Ny), cmap='binary', origin='lower', alpha=0.1)
            
            ax1.set_xticks(np.linspace(0,Nx,4))
            ax1.set_xticklabels([])
            ax1.set_xlabel(r'$x$', fontsize=12)
            ax1.set_yticks(np.linspace(0,Ny,4))
            ax1.set_yticklabels([])
            ax1.set_ylabel(r'$y$', fontsize=12)
            plt.colorbar(im1, location='right', shrink=0.97, aspect=20)

        
            ax2 = fig.add_subplot(gs[0,1])
            plt.title(r'$E_z^{scat}$')           
            im2 = plt.imshow(Ez_scat, vmin = -1, vmax = 1, cmap='seismic', origin='lower')
            plt.imshow(er, extent=(0, Nx, 0, Ny), cmap='binary', origin='lower', alpha=0.1)
            
            ax2.set_xticks(np.linspace(0,Nx,4))
            ax2.set_xticklabels([])
            ax2.set_xlabel(r'$x$', fontsize=12)
            ax2.set_yticks(np.linspace(0,Ny,4))
            ax2.set_yticklabels([])
##            ax2.set_ylabel(r'$y$', fontsize=12)

            plt.colorbar(im2, location='right', shrink=0.97, aspect=20)
       
##            plt.savefig(pictures+"/pic."+str(t).zfill(4)+".png")
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

















directory_info = 'data/information'
if not os.path.exists(directory_info):
    os.makedirs(directory_info)
file_name = "info_{}.txt".format(ratio)
file_path = os.path.join(directory_info, file_name)



'''computatio speed'''
# domain size Ny*Nx
# 6 fields (3 E fields and 3 H fields)
# 2 (scattered and total)
lattice_sites = 2 * 6 * Ny * Nx
time_steps = Time

mlups = lattice_sites * time_steps / ((t3 - t0) * 1e6)


import platform
import os
import subprocess

# Parse CPU information from /proc/cpuinfo
cpu_info = []
cache_size = []
with open("/proc/cpuinfo", "r") as file:
    for line in file:
        if line.strip():
            key, _, value = line.partition(":")
            if key.strip() == "model name":
                cpu_info.append(value.strip())
            if key.strip() == "cache size":
                cache_size.append(value.strip())


# Parse RAM information from /proc/meminfo
with open("/proc/meminfo", "r") as file:
    mem_info = {}
    for line in file:
        key, _, value = line.partition(":")
        mem_info[key.strip()] = value.strip()

# Get total and available memory in GB
total_memory = int(mem_info["MemTotal"].split()[0]) / 1024 / 1024  # Convert kB to GB
available_memory = int(mem_info["MemAvailable"].split()[0]) / 1024 / 1024  # Convert kB to GB



with open(file_path, "w") as file:
    file.write(f"Radius of the cylinder: {a}.\n")
    
    if (er2 > 1):
        file.write(f"Wavelength inside the scatterer: {wavelength * V2 / V1:.2f}.\n")
        
    file.write(f"Wavelength of the incident wave: {wavelength:.2f}.\n\n")
    
    file.write(f"Size of the computational domain: {Ny} * {Nx}.\n")
    file.write(f"Number of time steps: {Time}.\n\n")

    
    file.write(f"Number of parallel threads: {N}.\n")
    file.write(f"Total time taken: {t3 - t0:.2f} seconds.\n\n")
    file.write(f"Computation speed in MLUPS: {mlups:.2f}.\n\n")

    # CPU information
    file.write(f"CPU Model: {cpu_info[0]}\n")  # First processor
    file.write(f"Total logical cores: {len(cpu_info)}\n")
    file.write(f"Cache Size: {cache_size[0]}\n")  # Cache size of the first processor

    # system memory
    file.write(f"Total Memory: {total_memory:.2f} GB\n")
    file.write(f"Available Memory: {available_memory:.2f} GB\n\n")


    # python version
    file.write(f"Python Version: {platform.python_version()}\n")

    # Get GCC version
    try:
        gcc_version = subprocess.check_output(["gcc", "--version"], universal_newlines=True)
        gcc_version_line = gcc_version.splitlines()[0]
        file.write(f"GCC Version: {gcc_version_line}\n\n")
    except FileNotFoundError:
        file.write("GCC Compiler not found. Please install GCC to retrieve its version.\n\n")


    # Operating System
    file.write(f"Operating System: {platform.system()}\n")
    file.write(f"OS Version: {platform.version()}\n")
    file.write(f"OS Release: {platform.release()}\n\n")

    # Machine and Processor Info
    file.write(f"Machine: {platform.machine()}\n")
    file.write(f"Processor: {platform.processor()}\n")


print(f"Computation speed in MLUPS: {mlups:.2f}\n")

print('\n--------------------------------------------------------------------\n')

