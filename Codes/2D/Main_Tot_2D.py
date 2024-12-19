from scipy.integrate import simpson, trapezoid

from Module_Traction_2D import *
from Module_EM_Wave_2D import *
from Module_Parameters_2D import *
from Module_Shared_Lib_2D import *


t0 = time.time()


directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)





print("Number of time steps :", Time)




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

    
    if (t >= 0):
            
        '''source wave'''
        planeWaveTM(Ez, Hy, t, omega, xloc, ymin, ymax)

        '''collision and streaming (the 2 steps of LBM) when field is forced'''
        myclib.collForcingNode(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Ny, Nx, Q, xloc, ymin, ymax, N)
        myclib.streaming(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ny, Nx, Q, N)
    else:        
        '''collision and streaming (the 2 steps of LBM) when field is not forced'''
        myclib.collNotForcingNode(ex, ey, ez, hx, hy, hz, exb, eyb, ezb, hxb, hyb, hzb, Ex, Ey, Ez, Hx, Hy, Hz, er, mur, Ny, Nx, Q, N)
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

    FxIns[t] = (simpson(fxTop)   + simpson(fxRight)   + simpson(fxBottom)   + simpson(fxLeft))   / (a / ratio)
    FyIns[t] = (simpson(fyTop)   + simpson(fyRight)   + simpson(fyBottom)   + simpson(fyLeft))   / (a / ratio)
    TzIns[t] = (trapezoid(tzTop) + trapezoid(tzRight) + trapezoid(tzBottom) + trapezoid(tzLeft)) / (a / ratio)**2

    
    
    
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

Fx = (simpson(fxTop)   + simpson(fxRight)   + simpson(fxBottom)   + simpson(fxLeft))   / (a / ratio)
Fy = (simpson(fyTop)   + simpson(fyRight)   + simpson(fyBottom)   + simpson(fyLeft))   / (a / ratio)
Tz = (trapezoid(tzTop) + trapezoid(tzRight) + trapezoid(tzBottom) + trapezoid(tzLeft)) / (a / ratio)**2



t3 = time.time()
total_time = (t3 - t0) / 60
print(f"\nTotal time taken: {total_time:.2f} minutes\n")
###############################################################################################################

print(Fx)
print(Fy)
print(Tz)





np.save(directory+"/Ex_{}.txt".format(ratio), Ex_phasor)
np.save(directory+"/Ey_{}.txt".format(ratio), Ey_phasor)
np.save(directory+"/Ez_{}.txt".format(ratio), Ez_phasor)

np.save(directory+"/Hx_{}.txt".format(ratio), Hx_phasor)
np.save(directory+"/Hy_{}.txt".format(ratio), Hy_phasor)
np.save(directory+"/Hz_{}.txt".format(ratio), Hz_phasor)


#####################

np.save(directory+"/energy_{}.txt".format(ratio), U)

#####################

np.save(directory+"/FxIns_{}.txt".format(ratio), FxIns)
np.save(directory+"/FyIns_{}.txt".format(ratio), FyIns)
np.save(directory+"/TzIns_{}.txt".format(ratio), TzIns)




