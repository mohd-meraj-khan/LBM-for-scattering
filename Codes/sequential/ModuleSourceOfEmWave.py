import numpy as np



#################################################

def planeWaveTM(Ez, Hy, t, omega, xmin, xmax, ymin, ymax):
    Ez[0, ymin:ymax, xmin:xmax] = np.sin(omega * t)
    Hy[0, ymin:ymax, xmin:xmax] = - np.sin(omega * t)

#################################################




#################################################

def planeWaveTE(Ey, Hz, t, omega, xmin, xmax, ymin, ymax):
    Ey[0, ymin:ymax, xmin:xmax] = np.sin(omega * t)
    Hz[0, ymin:ymax, xmin:xmax] = np.sin(omega * t)

#################################################
