from Module_Traction_3D import *
from Module_EM_Wave_3D import *
from Module_Parameters_3D import *
from Module_Shared_Lib_3D import *


t0 = time.time()


directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)




print("Number of time steps :", int(Time))




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
Ex_phasor = np.zeros_like(Ex, dtype=complex)
Ey_phasor = np.zeros_like(Ey, dtype=complex)
Ez_phasor = np.zeros_like(Ez, dtype=complex)
Hx_phasor = np.zeros_like(Hx, dtype=complex)
Hy_phasor = np.zeros_like(Hy, dtype=complex)
Hz_phasor = np.zeros_like(Hz, dtype=complex)


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

        
    if (t >= 0):
            
        '''source wave'''
        planeWaveTM(EzI, HyI, t, omega, xloc, ymin, ymax, zmin, zmax)
        planeWaveTM(Ez, Hy, t, omega, xloc, ymin, ymax, zmin, zmax)

        '''calculation of scattered fields'''
        Ez_scat = Ez - EzI
        Hx_scat = Hx - HxI
        Hy_scat = Hy - HyI

##        Ez_scat[scatterer] = 0
##        Hx_scat[scatterer] = 0
##        Hy_scat[scatterer] = 0

        '''collision and streaming (the 2 steps of LBM) when field is forced'''
        myclib.collForcingNode(exI, eyI, ezI, hxI, hyI, hzI, exbI, eybI, ezbI, hxbI, hybI, hzbI, ExI, EyI, EzI, HxI, HyI, HzI, erI, murI, Nz, Ny, Nx, Q, xloc, ymin, ymax, zmin, zmax, N)
        myclib.collForcingNode(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Nz, Ny, Nx, Q, xloc, ymin, ymax, zmin, zmax, N)
        myclib.streaming(exI, eyI, ezI, hxI, hyI, hzI, exbI, eybI, ezbI, hxbI, hybI, hzbI, Nz, Ny, Nx, Q, N)
        myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Nz, Ny, Nx, Q, N)

    else:     
        '''collision and streaming (the 2 steps of LBM) when field is not forced'''
        myclib.collision(exI, eyI, ezI, hxI, hyI, hzI, exbI, eybI, ezbI, hxbI, hybI, hzbI, ExI, EyI, EzI, HxI, HyI, HzI, erI, murI, Nz, Ny, Nx, Q, N)
        myclib.collision(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Nz, Ny, Nx, Q, N)
        myclib.streaming(exI, eyI, ezI, hxI, hyI, hzI, exbI, eybI, ezbI, hxbI, hybI, hzbI, Nz, Ny, Nx, Q, N)
        myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Nz, Ny, Nx, Q, N)
         
    ###############################################################################################################
   
    U[t] = np.sum(0.5 * (er1 * (Ex**2 + Ey**2 + Ez**2) + mur1 * (Hx**2 + Hy**2 + Hz**2)))

    if (t >= Time - period):

        ###############################################################################################################
        
        '''converting from time domain to frequency domain'''
        Ex_phasor += Ex_scat * np.exp((0 - 1j) * omega * t) / period * 2
        Ey_phasor += Ey_scat * np.exp((0 - 1j) * omega * t) / period * 2
        Ez_phasor += Ez_scat * np.exp((0 - 1j) * omega * t) / period * 2

        Hx_phasor += Hx_scat * np.exp((0 - 1j) * omega * t) / period * 2
        Hy_phasor += Hy_scat * np.exp((0 - 1j) * omega * t) / period * 2
        Hz_phasor += Hz_scat * np.exp((0 - 1j) * omega * t) / period * 2
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


np.save(directory+"/ExScat_{}.npy".format(ratio), Ex_phasor)
np.save(directory+"/EyScat_{}.npy".format(ratio), Ey_phasor)
np.save(directory+"/EzScat_{}.npy".format(ratio), Ez_phasor)

np.save(directory+"/HxScat_{}.npy".format(ratio), Hx_phasor)
np.save(directory+"/HyScat_{}.npy".format(ratio), Hy_phasor)
np.save(directory+"/HzScat_{}.npy".format(ratio), Hz_phasor)

##########################################################

    
np.save(directory+"/energy_{}.npy".format(ratio), U)

