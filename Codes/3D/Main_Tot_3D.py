from scipy.integrate import simpson, trapezoid

from Module_Traction_3D import *
from Module_EM_Wave_3D import *
from Module_Parameters_3D import *
from Module_Shared_Lib_3D import *


t0 = time.time()


directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)





print("Number of time steps :", Time)




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
FzIns = np.zeros(Time)

TxIns = np.zeros(Time)
TyIns = np.zeros(Time)
TzIns = np.zeros(Time)


t1 = time.time()



for t in range(Time):
    
    

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
    '''total energy within the computational domain'''
    U[t] = np.sum(0.5 * (er1 * (Ex**2 + Ey**2 + Ez**2) + mur1 * (Hx**2 + Hy**2 + Hz**2)))



    '''instantaneous force and torque calculation'''
    fxTop    = fxInstantaneous(Ex_phasor[Top],    Ey_phasor[Top],    Ez_phasor[Top],    Hx_phasor[Top],    Hy_phasor[Top],    Hz_phasor[Top],    er1, mur1, nxTop,    nyTop,    nzTop)
    fxBottom = fxInstantaneous(Ex_phasor[Bottom], Ey_phasor[Bottom], Ez_phasor[Bottom], Hx_phasor[Bottom], Hy_phasor[Bottom], Hz_phasor[Bottom], er1, mur1, nxBottom, nyBottom, nzBottom)
    fxRight  = fxInstantaneous(Ex_phasor[Right],  Ey_phasor[Right],  Ez_phasor[Right],  Hx_phasor[Right],  Hy_phasor[Right],  Hz_phasor[Right],  er1, mur1, nxRight,  nyRight,  nzRight)
    fxLeft   = fxInstantaneous(Ex_phasor[Left],   Ey_phasor[Left],   Ez_phasor[Left],   Hx_phasor[Left],   Hy_phasor[Left],   Hz_phasor[Left],   er1, mur1, nxLeft,   nyLeft,   nzLeft)
    fxFront  = fxInstantaneous(Ex_phasor[Front],  Ey_phasor[Front],  Ez_phasor[Front],  Hx_phasor[Front],  Hy_phasor[Front],  Hz_phasor[Front],  er1, mur1, nxFront,  nyFront,  nzFront)
    fxBack   = fxInstantaneous(Ex_phasor[Back],   Ey_phasor[Back],   Ez_phasor[Back],   Hx_phasor[Back],   Hy_phasor[Back],   Hz_phasor[Back],   er1, mur1, nxBack,   nyBack,   nzBack)

    fyTop    = fyInstantaneous(Ex_phasor[Top],    Ey_phasor[Top],    Ez_phasor[Top],    Hx_phasor[Top],    Hy_phasor[Top],    Hz_phasor[Top],    er1, mur1, nxTop,    nyTop,    nzTop)
    fyBottom = fyInstantaneous(Ex_phasor[Bottom], Ey_phasor[Bottom], Ez_phasor[Bottom], Hx_phasor[Bottom], Hy_phasor[Bottom], Hz_phasor[Bottom], er1, mur1, nxBottom, nyBottom, nzBottom)
    fyRight  = fyInstantaneous(Ex_phasor[Right],  Ey_phasor[Right],  Ez_phasor[Right],  Hx_phasor[Right],  Hy_phasor[Right],  Hz_phasor[Right],  er1, mur1, nxRight,  nyRight,  nzRight)
    fyLeft   = fyInstantaneous(Ex_phasor[Left],   Ey_phasor[Left],   Ez_phasor[Left],   Hx_phasor[Left],   Hy_phasor[Left],   Hz_phasor[Left],   er1, mur1, nxLeft,   nyLeft,   nzLeft)
    fyFront  = fyInstantaneous(Ex_phasor[Front],  Ey_phasor[Front],  Ez_phasor[Front],  Hx_phasor[Front],  Hy_phasor[Front],  Hz_phasor[Front],  er1, mur1, nxFront,  nyFront,  nzFront)
    fyBack   = fyInstantaneous(Ex_phasor[Back],   Ey_phasor[Back],   Ez_phasor[Back],   Hx_phasor[Back],   Hy_phasor[Back],   Hz_phasor[Back],   er1, mur1, nxBack,   nyBack,   nzBack)

    fzTop    = fzInstantaneous(Ex_phasor[Top],    Ey_phasor[Top],    Ez_phasor[Top],    Hx_phasor[Top],    Hy_phasor[Top],    Hz_phasor[Top],    er1, mur1, nxTop,    nyTop,    nzTop)
    fzBottom = fzInstantaneous(Ex_phasor[Bottom], Ey_phasor[Bottom], Ez_phasor[Bottom], Hx_phasor[Bottom], Hy_phasor[Bottom], Hz_phasor[Bottom], er1, mur1, nxBottom, nyBottom, nzBottom)
    fzRight  = fzInstantaneous(Ex_phasor[Right],  Ey_phasor[Right],  Ez_phasor[Right],  Hx_phasor[Right],  Hy_phasor[Right],  Hz_phasor[Right],  er1, mur1, nxRight,  nyRight,  nzRight)
    fzLeft   = fzInstantaneous(Ex_phasor[Left],   Ey_phasor[Left],   Ez_phasor[Left],   Hx_phasor[Left],   Hy_phasor[Left],   Hz_phasor[Left],   er1, mur1, nxLeft,   nyLeft,   nzLeft)
    fzFront  = fzInstantaneous(Ex_phasor[Front],  Ey_phasor[Front],  Ez_phasor[Front],  Hx_phasor[Front],  Hy_phasor[Front],  Hz_phasor[Front],  er1, mur1, nxFront,  nyFront,  nzFront)
    fzBack   = fzInstantaneous(Ex_phasor[Back],   Ey_phasor[Back],   Ez_phasor[Back],   Hx_phasor[Back],   Hy_phasor[Back],   Hz_phasor[Back],   er1, mur1, nxBack,   nyBack,   nzBack)


    txTop    = (Top[1] - cy)    * fzTop    - (Top[0] - cz)    * fyTop
    txBottom = (Bottom[1] - cy) * fzBottom - (Bottom[0] - cz) * fyBottom
    txRight  = (Right[1] - cy)  * fzRight  - (Right[0] - cz)  * fyRight
    txLeft   = (Left[1] - cy)   * fzLeft   - (Left[0] - cz)   * fyLeft
    txFront  = (Front[1] - cy)  * fzFront  - (Front[0] - cz)  * fyFront
    txBack   = (Back[1] - cy)   * fzBack   - (Back[0] - cz)   * fyBack

    tyTop    = (Top[0] - cz)    * fxTop    - (Top[2] - cx)    * fzTop
    tyBottom = (Bottom[0] - cz) * fxBottom - (Bottom[2] - cx) * fzBottom
    tyRight  = (Right[0] - cz)  * fxRight  - (Right[2] - cx)  * fzRight
    tyLeft   = (Left[0] - cz)   * fxLeft   - (Left[2] - cx)   * fzLeft
    tyFront  = (Front[0] - cz)  * fxFront  - (Front[2] - cx)  * fzFront
    tyBack   = (Back[0] - cz)   * fxBack   - (Back[2] - cx)   * fzBack

    tzTop    = (Top[0] - cx)    * fyTop    - (Top[1] - cy)    * fxTop
    tzBottom = (Bottom[0] - cx) * fyBottom - (Bottom[1] - cy) * fxBottom
    tzRight  = (Right[0] - cx)  * fyRight  - (Right[1] - cy)  * fxRight
    tzLeft   = (Left[0] - cx)   * fyLeft   - (Left[1] - cy)   * fxLeft
    tzFront  = (Front[0] - cx)  * fyFront  - (Front[1] - cy)  * fxFront
    tzBack   = (Back[0] - cx)   * fyBack   - (Back[1] - cy)   * fxBack



    Fx = (simpson(fxTop) + simpson(fxBottom) + simpson(fxRight) + simpson(fxLeft) + simpson(fxFront) + simpson(fxBack)) / (a / ratio)
    Fy = (simpson(fyTop) + simpson(fyBottom) + simpson(fyRight) + simpson(fyLeft) + simpson(fyFront) + simpson(fyBack)) / (a / ratio)
    Fz = (simpson(fzTop) + simpson(fzBottom) + simpson(fzRight) + simpson(fzLeft) + simpson(fzFront) + simpson(fzBack)) / (a / ratio)

    Tx = (simpson(txTop) + simpson(txBottom) + simpson(txRight) + simpson(txLeft) + simpson(txFront) + simpson(txBack)) / (a / ratio)**2
    Ty = (simpson(tyTop) + simpson(tyBottom) + simpson(tyRight) + simpson(tyLeft) + simpson(tyFront) + simpson(tyBack)) / (a / ratio)**2
    Tz = (simpson(tzTop) + simpson(tzBottom) + simpson(tzRight) + simpson(tzLeft) + simpson(tzFront) + simpson(tzBack)) / (a / ratio)**2

    
    
    
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
fxBottom = fxAverage(Ex_phasor[Bottom], Ey_phasor[Bottom], Ez_phasor[Bottom], Hx_phasor[Bottom], Hy_phasor[Bottom], Hz_phasor[Bottom], er1, mur1, nxBottom, nyBottom, nzBottom)
fxRight  = fxAverage(Ex_phasor[Right],  Ey_phasor[Right],  Ez_phasor[Right],  Hx_phasor[Right],  Hy_phasor[Right],  Hz_phasor[Right],  er1, mur1, nxRight,  nyRight,  nzRight)
fxLeft   = fxAverage(Ex_phasor[Left],   Ey_phasor[Left],   Ez_phasor[Left],   Hx_phasor[Left],   Hy_phasor[Left],   Hz_phasor[Left],   er1, mur1, nxLeft,   nyLeft,   nzLeft)
fxFront  = fxAverage(Ex_phasor[Front],  Ey_phasor[Front],  Ez_phasor[Front],  Hx_phasor[Front],  Hy_phasor[Front],  Hz_phasor[Front],  er1, mur1, nxFront,  nyFront,  nzFront)
fxBack   = fxAverage(Ex_phasor[Back],   Ey_phasor[Back],   Ez_phasor[Back],   Hx_phasor[Back],   Hy_phasor[Back],   Hz_phasor[Back],   er1, mur1, nxBack,   nyBack,   nzBack)

fyTop    = fyAverage(Ex_phasor[Top],    Ey_phasor[Top],    Ez_phasor[Top],    Hx_phasor[Top],    Hy_phasor[Top],    Hz_phasor[Top],    er1, mur1, nxTop,    nyTop,    nzTop)
fyBottom = fyAverage(Ex_phasor[Bottom], Ey_phasor[Bottom], Ez_phasor[Bottom], Hx_phasor[Bottom], Hy_phasor[Bottom], Hz_phasor[Bottom], er1, mur1, nxBottom, nyBottom, nzBottom)
fyRight  = fyAverage(Ex_phasor[Right],  Ey_phasor[Right],  Ez_phasor[Right],  Hx_phasor[Right],  Hy_phasor[Right],  Hz_phasor[Right],  er1, mur1, nxRight,  nyRight,  nzRight)
fyLeft   = fyAverage(Ex_phasor[Left],   Ey_phasor[Left],   Ez_phasor[Left],   Hx_phasor[Left],   Hy_phasor[Left],   Hz_phasor[Left],   er1, mur1, nxLeft,   nyLeft,   nzLeft)
fyFront  = fyAverage(Ex_phasor[Front],  Ey_phasor[Front],  Ez_phasor[Front],  Hx_phasor[Front],  Hy_phasor[Front],  Hz_phasor[Front],  er1, mur1, nxFront,  nyFront,  nzFront)
fyBack   = fyAverage(Ex_phasor[Back],   Ey_phasor[Back],   Ez_phasor[Back],   Hx_phasor[Back],   Hy_phasor[Back],   Hz_phasor[Back],   er1, mur1, nxBack,   nyBack,   nzBack)

fzTop    = fzAverage(Ex_phasor[Top],    Ey_phasor[Top],    Ez_phasor[Top],    Hx_phasor[Top],    Hy_phasor[Top],    Hz_phasor[Top],    er1, mur1, nxTop,    nyTop,    nzTop)
fzBottom = fzAverage(Ex_phasor[Bottom], Ey_phasor[Bottom], Ez_phasor[Bottom], Hx_phasor[Bottom], Hy_phasor[Bottom], Hz_phasor[Bottom], er1, mur1, nxBottom, nyBottom, nzBottom)
fzRight  = fzAverage(Ex_phasor[Right],  Ey_phasor[Right],  Ez_phasor[Right],  Hx_phasor[Right],  Hy_phasor[Right],  Hz_phasor[Right],  er1, mur1, nxRight,  nyRight,  nzRight)
fzLeft   = fzAverage(Ex_phasor[Left],   Ey_phasor[Left],   Ez_phasor[Left],   Hx_phasor[Left],   Hy_phasor[Left],   Hz_phasor[Left],   er1, mur1, nxLeft,   nyLeft,   nzLeft)
fzFront  = fzAverage(Ex_phasor[Front],  Ey_phasor[Front],  Ez_phasor[Front],  Hx_phasor[Front],  Hy_phasor[Front],  Hz_phasor[Front],  er1, mur1, nxFront,  nyFront,  nzFront)
fzBack   = fzAverage(Ex_phasor[Back],   Ey_phasor[Back],   Ez_phasor[Back],   Hx_phasor[Back],   Hy_phasor[Back],   Hz_phasor[Back],   er1, mur1, nxBack,   nyBack,   nzBack)


txTop    = (Top[1] - cy)    * fzTop    - (Top[0] - cz)    * fyTop
txBottom = (Bottom[1] - cy) * fzBottom - (Bottom[0] - cz) * fyBottom
txRight  = (Right[1] - cy)  * fzRight  - (Right[0] - cz)  * fyRight
txLeft   = (Left[1] - cy)   * fzLeft   - (Left[0] - cz)   * fyLeft
txFront  = (Front[1] - cy)  * fzFront  - (Front[0] - cz)  * fyFront
txBack   = (Back[1] - cy)   * fzBack   - (Back[0] - cz)   * fyBack

tyTop    = (Top[0] - cz)    * fxTop    - (Top[2] - cx)    * fzTop
tyBottom = (Bottom[0] - cz) * fxBottom - (Bottom[2] - cx) * fzBottom
tyRight  = (Right[0] - cz)  * fxRight  - (Right[2] - cx)  * fzRight
tyLeft   = (Left[0] - cz)   * fxLeft   - (Left[2] - cx)   * fzLeft
tyFront  = (Front[0] - cz)  * fxFront  - (Front[2] - cx)  * fzFront
tyBack   = (Back[0] - cz)   * fxBack   - (Back[2] - cx)   * fzBack

tzTop    = (Top[0] - cx)    * fyTop    - (Top[1] - cy)    * fxTop
tzBottom = (Bottom[0] - cx) * fyBottom - (Bottom[1] - cy) * fxBottom
tzRight  = (Right[0] - cx)  * fyRight  - (Right[1] - cy)  * fxRight
tzLeft   = (Left[0] - cx)   * fyLeft   - (Left[1] - cy)   * fxLeft
tzFront  = (Front[0] - cx)  * fyFront  - (Front[1] - cy)  * fxFront
tzBack   = (Back[0] - cx)   * fyBack   - (Back[1] - cy)   * fxBack



Fx = (simpson(fxTop) + simpson(fxBottom) + simpson(fxRight) + simpson(fxLeft) + simpson(fxFront) + simpson(fxBack)) / (a / ratio)
Fy = (simpson(fyTop) + simpson(fyBottom) + simpson(fyRight) + simpson(fyLeft) + simpson(fyFront) + simpson(fyBack)) / (a / ratio)
Fz = (simpson(fzTop) + simpson(fzBottom) + simpson(fzRight) + simpson(fzLeft) + simpson(fzFront) + simpson(fzBack)) / (a / ratio)

Tx = (simpson(txTop) + simpson(txBottom) + simpson(txRight) + simpson(txLeft) + simpson(txFront) + simpson(txBack)) / (a / ratio)**2
Ty = (simpson(tyTop) + simpson(tyBottom) + simpson(tyRight) + simpson(tyLeft) + simpson(tyFront) + simpson(tyBack)) / (a / ratio)**2
Tz = (simpson(tzTop) + simpson(tzBottom) + simpson(tzRight) + simpson(tzLeft) + simpson(tzFront) + simpson(tzBack)) / (a / ratio)**2



t3 = time.time()
total_time = (t3 - t0) / 60
print(f"\nTotal time taken: {total_time:.2f} minutes\n")
###############################################################################################################

print(Fx)
print(Fy)
print(Fz)

print(Tx)
print(Ty)
print(Tz)



np.save(directory+"/Ex_{}.npy".format(ratio), Ex_phasor)
np.save(directory+"/Ey_{}.npy".format(ratio), Ey_phasor)
np.save(directory+"/Ez_{}.npy".format(ratio), Ez_phasor)

np.save(directory+"/Hx_{}.npy".format(ratio), Hx_phasor)
np.save(directory+"/Hy_{}.npy".format(ratio), Hy_phasor)
np.save(directory+"/Hz_{}.npy".format(ratio), Hz_phasor)



##########################################################

np.save(directory+"/energy_{}.npy".format(ratio), U)


##########################################################

np.save(directory+"/FxIns_{}.npy".format(ratio), FxIns)
np.save(directory+"/FyIns_{}.npy".format(ratio), FyIns)
np.save(directory+"/FzIns_{}.npy".format(ratio), FzIns)

np.save(directory+"/TxIns_{}.npy".format(ratio), TxIns)
np.save(directory+"/TyIns_{}.npy".format(ratio), TyIns)
np.save(directory+"/TzIns_{}.npy".format(ratio), TzIns)




