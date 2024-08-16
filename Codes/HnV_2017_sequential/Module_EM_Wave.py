import numpy as np



#################################################

def planeWaveTM(Ez, Hy, t, omega, x, ymin, ymax):
    Ez[ymin:ymax, x] = np.sin(omega * t)
    Hy[ymin:ymax, x] = - np.sin(omega * t)

#################################################




#################################################

def planeWaveTE(Ey, Hz, t, omega, x, ymin, ymax):
    Ey[ymin:ymax, x] = np.sin(omega * t)
    Hz[ymin:ymax, x] = np.sin(omega * t)

#################################################
