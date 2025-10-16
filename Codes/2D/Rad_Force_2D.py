import numpy as np

from scipy.integrate import simpson

import sys
import os

from Module_Traction_2D import *
from Module_Parameters_2D import *





directory = 'data/total_field'


force = 'data'
if not os.path.exists(force):
    os.makedirs(force)




Ex = np.load(directory+'/ExTot_{}_{}.npy'.format(er2, ratio))
Ey = np.load(directory+'/EyTot_{}_{}.npy'.format(er2, ratio))
Ez = np.load(directory+'/EzTot_{}_{}.npy'.format(er2, ratio))

Hx = np.load(directory+'/HxTot_{}_{}.npy'.format(er2, ratio))
Hy = np.load(directory+'/HyTot_{}_{}.npy'.format(er2, ratio))
Hz = np.load(directory+'/HzTot_{}_{}.npy'.format(er2, ratio))




'''center of the scatterer'''
cx = Nx//2 + 0.5
cy = Ny//2 + 0.5

'''half-width of the bounding box'''
w = int(np.round(1.5*a))
'''box surrounding the scatterer'''
Top    = int(cy + w), np.arange(int(cx - w), int(cx + w))
Right  = np.arange(int(cy - w), int(cy + w)), int(cx + w)
Bottom = int(cy - w), np.arange(int(cx - w), int(cx + w))
Left   = np.arange(int(cy - w), int(cy + w)), int(cx - w)

'''unit normal vectors at the perimeter of the bounding box'''
nxTop, nxRight, nxBottom, nxLeft = 0, 1, 0, -1
nyTop, nyRight, nyBottom, nyLeft = 1, 0, -1, 0
nzTop, nzRight, nzBottom, nzLeft = 0, 0,  0, 0


'''force and torque calculation averaged over one time period'''
fxTop    = fxAverage(Ex[Top], Ey[Top], Ez[Top], Hx[Top], Hy[Top], Hz[Top], er1, mur1, nxTop, nyTop, nzTop)
fxRight  = fxAverage(Ex[Right], Ey[Right], Ez[Right], Hx[Right], Hy[Right], Hz[Right], er1, mur1, nxRight, nyRight, nzRight)
fxBottom = fxAverage(Ex[Bottom], Ey[Bottom], Ez[Bottom], Hx[Bottom], Hy[Bottom], Hz[Bottom], er1, mur1, nxBottom, nyBottom, nzBottom)
fxLeft   = fxAverage(Ex[Left], Ey[Left], Ez[Left], Hx[Left], Hy[Left], Hz[Left], er1, mur1, nxLeft, nyLeft, nzLeft)

fyTop    = fyAverage(Ex[Top],    Ey[Top],    Ez[Top],    Hx[Top],    Hy[Top],    Hz[Top],    er1, mur1, nxTop,    nyTop,    nzTop)
fyRight  = fyAverage(Ex[Right],  Ey[Right],  Ez[Right],  Hx[Right],  Hy[Right],  Hz[Right],  er1, mur1, nxRight,  nyRight,  nzRight)
fyBottom = fyAverage(Ex[Bottom], Ey[Bottom], Ez[Bottom], Hx[Bottom], Hy[Bottom], Hz[Bottom], er1, mur1, nxBottom, nyBottom, nzBottom)
fyLeft   = fyAverage(Ex[Left],   Ey[Left],   Ez[Left],   Hx[Left],   Hy[Left],   Hz[Left],   er1, mur1, nxLeft,   nyLeft,   nzLeft)

tzTop    = - fxTop    * (Top[0] - cy)    + fyTop    * (Top[1] - cx)
tzRight  = - fxRight  * (Right[0] - cx)  + fyRight  * (Right[1] - cy)
tzBottom = - fxBottom * (Bottom[0] - cy) + fyBottom * (Bottom[1] - cx)
tzLeft   = - fxLeft   * (Left[0] - cx)   + fyLeft   * (Left[1] - cy)


Fx = (simpson(fxTop) + simpson(fxRight) + simpson(fxBottom) + simpson(fxLeft)) / (a / ratio)
Fy = (simpson(fyTop) + simpson(fyRight) + simpson(fyBottom) + simpson(fyLeft)) / (a / ratio)
Tz = (simpson(tzTop) + simpson(tzRight) + simpson(tzBottom) + simpson(tzLeft)) / (a / ratio)**2



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




