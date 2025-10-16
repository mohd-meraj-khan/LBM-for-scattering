from scipy.integrate import simps

from Module_Traction_2D import *
from Module_EM_Wave_2D import *
from Module_Parameters_2D import *
from Module_Shared_Lib_2D import *


t0 = time.time()


field = 'data/total_field'
if not os.path.exists(field):
    os.makedirs(field)

force = 'data'
if not os.path.exists(force):
    os.makedirs(force)

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
def initialize_field(Ny=10, Nx=10):
    return np.zeros((Ny, Nx), dtype=np.float32, order='C')

Ex, Ey, Ez, Hx, Hy, Hz = [initialize_field(Ny, Nx) for _ in range(6)]


'''initializing the distribution functions of electric and magnetic fields'''
def initilize_dis_func(Ny=10, Nx=10, Q=7):
    return np.zeros((Ny, Nx, Q), dtype=np.float32, order='C')

ex, ey, ez, hx, hy, hz = [initilize_dis_func(Ny, Nx, Q) for _ in range(6)]
exb, eyb, ezb, hxb, hyb, hzb = [initilize_dis_func(Ny, Nx, Q) for _ in range(6)]

###############################################################################################################



'''initilizing the variables for frequency domain fields'''
Ex_phasor = np.zeros_like(Ex, dtype=complex)
Ey_phasor = np.zeros_like(Ey, dtype=complex)
Ez_phasor = np.zeros_like(Ez, dtype=complex)
Hx_phasor = np.zeros_like(Hx, dtype=complex)
Hy_phasor = np.zeros_like(Hy, dtype=complex)
Hz_phasor = np.zeros_like(Hz, dtype=complex)


U = np.zeros(Time)

FxIns = np.zeros(Time)
FyIns = np.zeros(Time)
TzIns = np.zeros(Time)


t1 = time.time()



for t in range(Time):
    
    

    #################################################################################################################
    ########                                         LBM CALCULATION                                          #######
    #################################################################################################################

    '''initialization of macroscopic fields'''
    myclib.initializeField(Ex, Ey, Ez, Hx, Hy, Hz, Ny, Nx, N)


    '''computation of macroscopic fields from particle populations'''
    myclib.macroField(ex, er, Ex, Ny, Nx, Q, N)
    myclib.macroField(ey, er, Ey, Ny, Nx, Q, N)
    myclib.macroField(ez, er, Ez, Ny, Nx, Q, N)

    myclib.macroField(hx, mur, Hx, Ny, Nx, Q, N)
    myclib.macroField(hy, mur, Hy, Ny, Nx, Q, N)
    myclib.macroField(hz, mur, Hz, Ny, Nx, Q, N)

            
    '''source wave'''
    planeWaveTM(Ez, Hy, t, omega, xloc, ymin, ymax)


    '''electric field is zero inside PEC'''
##    Ex[scatterer] = 0
##    Ey[scatterer] = 0
##    Ez[scatterer] = 0


    '''collision and streaming (the 2 steps of LBM) when field is forced'''
    myclib.collForcingNode(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Ny, Nx, Q, xloc, ymin, ymax, N)
    myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ny, Nx, Q, N)


    ###############################################################################################################
    '''total energy within the computational domain'''
    U[t] = np.sum(0.5 * (er1 * (Ex**2 + Ey**2 + Ez**2) + mur1 * (Hx**2 + Hy**2 + Hz**2)))



    '''instantaneous force and torque calculation'''
    fxTop    = fxInstantaneous(Ex[Top],    Ey[Top],    Ez[Top],    Hx[Top],    Hy[Top],    Hz[Top],    er1, mur1, nxTop,    nyTop,    nzTop)
    fxRight  = fxInstantaneous(Ex[Right],  Ey[Right],  Ez[Right],  Hx[Right],  Hy[Right],  Hz[Right],  er1, mur1, nxRight,  nyRight,  nzRight)
    fxBottom = fxInstantaneous(Ex[Bottom], Ey[Bottom], Ez[Bottom], Hx[Bottom], Hy[Bottom], Hz[Bottom], er1, mur1, nxBottom, nyBottom, nzBottom)
    fxLeft   = fxInstantaneous(Ex[Left],   Ey[Left],   Ez[Left],   Hx[Left],   Hy[Left],   Hz[Left],   er1, mur1, nxLeft,   nyLeft,   nzLeft)

    fyTop    = fyInstantaneous(Ex[Top],    Ey[Top],    Ez[Top],    Hx[Top],    Hy[Top],    Hz[Top],    er1, mur1, nxTop,    nyTop,    nzTop)
    fyRight  = fyInstantaneous(Ex[Right],  Ey[Right],  Ez[Right],  Hx[Right],  Hy[Right],  Hz[Right],  er1, mur1, nxRight,  nyRight,  nzRight)
    fyBottom = fyInstantaneous(Ex[Bottom], Ey[Bottom], Ez[Bottom], Hx[Bottom], Hy[Bottom], Hz[Bottom], er1, mur1, nxBottom, nyBottom, nzBottom)
    fyLeft   = fyInstantaneous(Ex[Left],   Ey[Left],   Ez[Left],   Hx[Left],   Hy[Left],   Hz[Left],   er1, mur1, nxLeft,   nyLeft,   nzLeft)

    tzTop    = - fxTop    * (Top[0] - cy)    + fyTop    * (Top[1] - cx)
    tzRight  = - fxRight  * (Right[0] - cx)  + fyRight  * (Right[1] - cy)
    tzBottom = - fxBottom * (Bottom[0] - cy) + fyBottom * (Bottom[1] - cx)
    tzLeft   = - fxLeft   * (Left[0] - cx)   + fyLeft   * (Left[1] - cy)

    FxIns[t] = (simps(fxTop) + simps(fxRight) + simps(fxBottom) + simps(fxLeft)) / (a / ratio)
    FyIns[t] = (simps(fyTop) + simps(fyRight) + simps(fyBottom) + simps(fyLeft)) / (a / ratio)
    TzIns[t] = (simps(tzTop) + simps(tzRight) + simps(tzBottom) + simps(tzLeft)) / (a / ratio)**2

    
    
    
    if (t >= Time - period):
        
        '''converting from time domain to frequency domain'''
        Ex_phasor += Ex * np.exp((0 - 1j) * omega * t) / period * 2
        Ey_phasor += Ey * np.exp((0 - 1j) * omega * t) / period * 2
        Ez_phasor += Ez * np.exp((0 - 1j) * omega * t) / period * 2

        Hx_phasor += Hx * np.exp((0 - 1j) * omega * t) / period * 2
        Hy_phasor += Hy * np.exp((0 - 1j) * omega * t) / period * 2
        Hz_phasor += Hz * np.exp((0 - 1j) * omega * t) / period * 2
        ###############################################################################################################


###############################################################################################################   
    t2 = time.time()
        
    if (t > 0 and t%100 == 0):
        remaining_time = (t2 - t1) * (int(Time) - t) / (t*60)
        print(f"Approximate time left: {remaining_time:.2f} minutes", end="\r")


'''force and torque calculation averaged over one time period'''
fxTop    = fxAverage(Ex_phasor[Top],    Ey_phasor[Top],    Ez_phasor[Top],    Hx_phasor[Top],    Hy_phasor[Top],    Hz_phasor[Top],    er1, mur1, nxTop,    nyTop,    nzTop)
fxRight  = fxAverage(Ex_phasor[Right],  Ey_phasor[Right],  Ez_phasor[Right],  Hx_phasor[Right],  Hy_phasor[Right],  Hz_phasor[Right],  er1, mur1, nxRight,  nyRight,  nzRight)
fxBottom = fxAverage(Ex_phasor[Bottom], Ey_phasor[Bottom], Ez_phasor[Bottom], Hx_phasor[Bottom], Hy_phasor[Bottom], Hz_phasor[Bottom], er1, mur1, nxBottom, nyBottom, nzBottom)
fxLeft   = fxAverage(Ex_phasor[Left],   Ey_phasor[Left],   Ez_phasor[Left],   Hx_phasor[Left],   Hy_phasor[Left],   Hz_phasor[Left],   er1, mur1, nxLeft,   nyLeft,   nzLeft)

fyTop    = fyAverage(Ex_phasor[Top],    Ey_phasor[Top],    Ez_phasor[Top],    Hx_phasor[Top],    Hy_phasor[Top],    Hz_phasor[Top],    er1, mur1, nxTop,    nyTop,    nzTop)
fyRight  = fyAverage(Ex_phasor[Right],  Ey_phasor[Right],  Ez_phasor[Right],  Hx_phasor[Right],  Hy_phasor[Right],  Hz_phasor[Right],  er1, mur1, nxRight,  nyRight,  nzRight)
fyBottom = fyAverage(Ex_phasor[Bottom], Ey_phasor[Bottom], Ez_phasor[Bottom], Hx_phasor[Bottom], Hy_phasor[Bottom], Hz_phasor[Bottom], er1, mur1, nxBottom, nyBottom, nzBottom)
fyLeft   = fyAverage(Ex_phasor[Left],   Ey_phasor[Left],   Ez_phasor[Left],   Hx_phasor[Left],   Hy_phasor[Left],   Hz_phasor[Left],   er1, mur1, nxLeft,   nyLeft,   nzLeft)

tzTop    = - fxTop    * (Top[0] - cy)    + fyTop    * (Top[1] - cx)
tzRight  = - fxRight  * (Right[0] - cx)  + fyRight  * (Right[1] - cy)
tzBottom = - fxBottom * (Bottom[0] - cy) + fyBottom * (Bottom[1] - cx)
tzLeft   = - fxLeft   * (Left[0] - cx)   + fyLeft   * (Left[1] - cy)


Fx = (simps(fxTop) + simps(fxRight) + simps(fxBottom) + simps(fxLeft)) / (a / ratio)
Fy = (simps(fyTop) + simps(fyRight) + simps(fyBottom) + simps(fyLeft)) / (a / ratio)
Tz = (simps(tzTop) + simps(tzRight) + simps(tzBottom) + simps(tzLeft)) / (a / ratio)**2



t3 = time.time()
total_time = (t3 - t0) / 60
print(f"\nTotal time taken: {total_time:.2f} minutes\n")
###############################################################################################################

print(f'Fx: {Fx}')
print(f'Fy: {Fy}')
print(f'Tz: {Tz}')

Fx_avg = []
Fy_avg = []
Tz_avg = []

Fx_avg.append(Fx)
Fy_avg.append(Fy)
Tz_avg.append(Tz)

fx = open(force+"/Fx_er_{}.txt".format(er2), "a")
np.savetxt(fx, Fx_avg)
fx.close()

fy = open(force+"/Fy_er_{}.txt".format(er2), "a")
np.savetxt(fy, Fy_avg)
fy.close()

fz = open(force+"/Tz_er_{}.txt".format(er2), "a")
np.savetxt(fz, Tz_avg)
fz.close()



np.save(field+"/ExTot_{}_{}.npy".format(er2, ratio), Ex_phasor)
np.save(field+"/EyTot_{}_{}.npy".format(er2, ratio), Ey_phasor)
np.save(field+"/EzTot_{}_{}.npy".format(er2, ratio), Ez_phasor)

np.save(field+"/HxTot_{}_{}.npy".format(er2, ratio), Hx_phasor)
np.save(field+"/HyTot_{}_{}.npy".format(er2, ratio), Hy_phasor)
np.save(field+"/HzTot_{}_{}.npy".format(er2, ratio), Hz_phasor)


#####################

np.save(MA+"/energy_{}_{}.npy".format(er2, ratio), U)

#####################

np.save(MA+"/FxIns_{}_{}.npy".format(er2, ratio), FxIns)
np.save(MA+"/FyIns_{}_{}.npy".format(er2, ratio), FyIns)
np.save(MA+"/TzIns_{}_{}.npy".format(er2, ratio), TzIns)



###################################################



directory_info = 'data/information'
if not os.path.exists(directory_info):
    os.makedirs(directory_info)

    
file_name = "info_{}_{}.txt".format(er2, ratio)
file_path = os.path.join(directory_info, file_name)



'''computatio speed'''
# domain size Ny*Nx
# 6 fields (3 E fields and 3 H fields)
lattice_sites = 6 * Ny * Nx
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


