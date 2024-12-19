import numpy as np



#################################################

def planeWaveTM(Ez, Hy, t, omega, x, ymin, ymax):
    Ez[0, ymin:ymax, x] = np.sin(omega * t)
    Hy[0, ymin:ymax, x] = - np.sin(omega * t)

#################################################




#################################################

def planeWaveTE(Ey, Hz, t, omega, x, ymin, ymax):
    Ey[0, ymin:ymax, x] = np.sin(omega * t)
    Hz[0, ymin:ymax, x] = np.sin(omega * t)

#################################################
