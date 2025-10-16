import numpy as np




'''TMz plane wave'''
def planeWaveTM(Ez, Hy, t, omega, x, ymin, ymax, zmin, zmax):
    Ez[zmin:zmax, ymin:ymax, x] =  np.sin(omega * t)
    Hy[zmin:zmax, ymin:ymax, x] = - np.sin(omega * t)


'''TEz plane wave'''
def planeWaveTE(Ey, Hz, t, omega, x, ymin, ymax, zmin, zmax):
    Ey[zmin:zmax, ymin:ymax, x] = np.sin(omega * t)
    Hz[zmin:zmax, ymin:ymax, x] = np.sin(omega * t)


