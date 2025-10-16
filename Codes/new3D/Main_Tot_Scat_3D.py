from Module_Traction_3D import *
from Module_EM_Wave_3D import *
from Module_Parameters_3D import *
from Module_Shared_Lib_3D import *


t0 = time.time()


directory_scattered = 'data/scattered_field'
if not os.path.exists(directory_scattered):
    os.makedirs(directory_scattered)

directory_total = 'data/total_field'
if not os.path.exists(directory_total):
    os.makedirs(directory_total)


MA = 'data/moving_average'
if not os.path.exists(MA):
    os.makedirs(MA)


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
def initialize_field(Nz=10, Ny=10, Nx=10):
    return np.zeros((Nz, Ny, Nx), dtype=np.float32, order='C')

ExI, EyI, EzI, HxI, HyI, HzI = [initialize_field(Nz, Ny, Nx) for _ in range(6)]
Ex, Ey, Ez, Hx, Hy, Hz = [initialize_field(Nz, Ny, Nx) for _ in range(6)]
Ex_scat, Ey_scat, Ez_scat, Hx_scat, Hy_scat, Hz_scat = [initialize_field(Nz, Ny, Nx) for _ in range(6)]


'''initializing the distribution functions of electric and magnetic fields'''
def initilize_dis_func(Nz=10, Ny=10, Nx=10, Q=7):
    return np.zeros((Nz, Ny, Nx, Q), dtype=np.float32, order='C')

exI, eyI, ezI, hxI, hyI, hzI = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]
exbI, eybI, ezbI, hxbI, hybI, hzbI = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]
ex, ey, ez, hx, hy, hz = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]
exb, eyb, ezb, hxb, hyb, hzb = [initilize_dis_func(Nz, Ny, Nx, Q) for _ in range(6)]

'''initializing the material properties for incident fields'''
erI, murI = initialize_material_properties(er1, mur1, Nz, Ny, Nx)
###############################################################################################################



'''initilizing the variables for frequency domain fields'''
ExScat = np.zeros_like(Ex, dtype=complex)
EyScat = np.zeros_like(Ey, dtype=complex)
EzScat = np.zeros_like(Ez, dtype=complex)
HxScat = np.zeros_like(Hx, dtype=complex)
HyScat = np.zeros_like(Hy, dtype=complex)
HzScat = np.zeros_like(Hz, dtype=complex)

ExTot  = np.zeros_like(Ex, dtype=complex)
EyTot  = np.zeros_like(Ey, dtype=complex)
EzTot  = np.zeros_like(Ez, dtype=complex)
HxTot  = np.zeros_like(Hx, dtype=complex)
HyTot  = np.zeros_like(Hy, dtype=complex)
HzTot  = np.zeros_like(Hz, dtype=complex)


U = np.zeros(Time)

t1 = time.time()





for t in range(int(Time)):



    #################################################################################################################
    ########                                         LBM CALCULATION                                          #######
    #################################################################################################################

    '''initialization of macroscopic fields'''
    myclib.initializeField(ExI, EyI, EzI, HxI, HyI, HzI, Nz, Ny, Nx, N)
    myclib.initializeField(Ex, Ey, Ez, Hx, Hy, Hz, Nz, Ny, Nx, N)

    '''computation of macroscopic fields from distribution function'''
    myclib.macroField(exI, erI, ExI, Nz, Ny, Nx, Q, N)
    myclib.macroField(eyI, erI, EyI, Nz, Ny, Nx, Q, N)
    myclib.macroField(ezI, erI, EzI, Nz, Ny, Nx, Q, N)

    myclib.macroField(hxI, murI, HxI, Nz, Ny, Nx, Q, N)
    myclib.macroField(hyI, murI, HyI, Nz, Ny, Nx, Q, N)
    myclib.macroField(hzI, murI, HzI, Nz, Ny, Nx, Q, N)


    myclib.macroField(ex, er, Ex, Nz, Ny, Nx, Q, N)
    myclib.macroField(ey, er, Ey, Nz, Ny, Nx, Q, N)
    myclib.macroField(ez, er, Ez, Nz, Ny, Nx, Q, N)

    myclib.macroField(hx, mur, Hx, Nz, Ny, Nx, Q, N)
    myclib.macroField(hy, mur, Hy, Nz, Ny, Nx, Q, N)
    myclib.macroField(hz, mur, Hz, Nz, Ny, Nx, Q, N)

        
            
    '''source wave'''
    planeWaveTM(EzI, HyI, t, omega, xloc, ymin, ymax, zmin, zmax)
    planeWaveTM(Ez, Hy, t, omega, xloc, ymin, ymax, zmin, zmax)

    '''calculation of scattered fields'''
    Ex_scat = Ex - ExI
    Ey_scat = Ey - EyI
    Ez_scat = Ez - EzI
    Hx_scat = Hx - HxI
    Hy_scat = Hy - HyI
    Hz_scat = Hz - HzI
    

    '''electric field is zero inside PEC'''
##    Ex[scatterer] = 0
##    Ey[scatterer] = 0
##    Ez[scatterer] = 0


    '''collision and streaming (the 2 steps of LBM) when field is forced'''
    myclib.collForcingNode(exI, eyI, ezI, hxI, hyI, hzI, exbI, eybI, ezbI, hxbI, hybI, hzbI, ExI, EyI, EzI, HxI, HyI, HzI, erI, murI, Nz, Ny, Nx, Q, xloc, ymin, ymax, zmin, zmax, N)
    myclib.collForcingNode(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Nz, Ny, Nx, Q, xloc, ymin, ymax, zmin, zmax, N)
    myclib.streaming(exI, eyI, ezI, hxI, hyI, hzI, exbI, eybI, ezbI, hxbI, hybI, hzbI, Nz, Ny, Nx, Q, N)
    myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Nz, Ny, Nx, Q, N)

            
    ###############################################################################################################
   
    U[t] = np.sum(0.5 * (er1 * (Ex**2 + Ey**2 + Ez**2) + mur1 * (Hx**2 + Hy**2 + Hz**2)))

    if (t >= Time - period):

        ###############################################################################################################
        
        '''converting from time domain to frequency domain'''
        ExScat += Ex_scat * np.exp((0 - 1j) * omega * t) / period * 2
        EyScat += Ey_scat * np.exp((0 - 1j) * omega * t) / period * 2
        EzScat += Ez_scat * np.exp((0 - 1j) * omega * t) / period * 2
        HxScat += Hx_scat * np.exp((0 - 1j) * omega * t) / period * 2
        HyScat += Hy_scat * np.exp((0 - 1j) * omega * t) / period * 2
        HzScat += Hz_scat * np.exp((0 - 1j) * omega * t) / period * 2

        ExTot += Ex * np.exp((0 - 1j) * omega * t) / period * 2
        EyTot += Ey * np.exp((0 - 1j) * omega * t) / period * 2
        EzTot += Ez * np.exp((0 - 1j) * omega * t) / period * 2
        HxTot += Hx * np.exp((0 - 1j) * omega * t) / period * 2
        HyTot += Hy * np.exp((0 - 1j) * omega * t) / period * 2
        HzTot += Hz * np.exp((0 - 1j) * omega * t) / period * 2
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


np.save(directory_scattered+"/ExScat_{}_{}.npy".format(er2, ratio), ExScat)
np.save(directory_scattered+"/EyScat_{}_{}.npy".format(er2, ratio), EyScat)
np.save(directory_scattered+"/EzScat_{}_{}.npy".format(er2, ratio), EzScat)
np.save(directory_scattered+"/HxScat_{}_{}.npy".format(er2, ratio), HxScat)
np.save(directory_scattered+"/HyScat_{}_{}.npy".format(er2, ratio), HyScat)
np.save(directory_scattered+"/HzScat_{}_{}.npy".format(er2, ratio), HzScat)

np.save(directory_total+"/ExTot_{}_{}.npy".format(er2, ratio), ExTot)
np.save(directory_total+"/EyTot_{}_{}.npy".format(er2, ratio), EyTot)
np.save(directory_total+"/EzTot_{}_{}.npy".format(er2, ratio), EzTot)
np.save(directory_total+"/HxTot_{}_{}.npy".format(er2, ratio), HxTot)
np.save(directory_total+"/HyTot_{}_{}.npy".format(er2, ratio), HyTot)
np.save(directory_total+"/HzTot_{}_{}.npy".format(er2, ratio), HzTot)


#####################

np.save(MA+"/energy_{}_{}.npy".format(er2, ratio), U)

#####################



directory_info = 'data/information'
if not os.path.exists(directory_info):
    os.makedirs(directory_info)
file_name = "info_{}_{}.txt".format(er2, ratio)
file_path = os.path.join(directory_info, file_name)



'''computatio speed'''
# domain size Ny*Nx
# 6 fields (3 E fields and 3 H fields)
# 2 (scattered and total)
lattice_sites = 2 * 6 * Nz * Ny * Nx
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


