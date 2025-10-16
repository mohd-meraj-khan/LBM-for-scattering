import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from scipy.integrate import simpson

import sys
import os

from Module_Traction_2D import *
from Module_Parameters_2D import *





directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)

plots = 'plots'
if not os.path.exists(plots):
    os.makedirs(plots)




Ex = np.genfromtxt(directory+'/Ex_{}.txt'.format(ratio), dtype=complex)
Ey = np.genfromtxt(directory+'/Ey_{}.txt'.format(ratio), dtype=complex)
Ez = np.genfromtxt(directory+'/Ez_{}.txt'.format(ratio), dtype=complex)

Hx = np.genfromtxt(directory+'/Hx_{}.txt'.format(ratio), dtype=complex)
Hy = np.genfromtxt(directory+'/Hy_{}.txt'.format(ratio), dtype=complex)
Hz = np.genfromtxt(directory+'/Hz_{}.txt'.format(ratio), dtype=complex)




'''center of the scatterer'''
cx = Nx//2 + 0.5
cy = Ny//2 + 0.5
cz = Nz//2 + 0.5

'''half-width of the square bounding box'''
w = int(np.round(1.5*a))

'''box surrounding the scatterer'''
Top    = np.arange(int(cz - w), int(cz + w)), int(cy + w), np.arange(int(cx - w), int(cx + w))
Bottom = np.arange(int(cz - w), int(cz + w)), int(cy - w), np.arange(int(cx - w), int(cx + w))
Right  = np.arange(int(cz - w), int(cz + w)), np.arange(int(cy - w), int(cy + w)), int(cx + w)
Left   = np.arange(int(cz - w), int(cz + w)), np.arange(int(cy - w), int(cy + w)), int(cx - w)
Front  = int(cz + w), np.arange(int(cy - w), int(cy + w)), np.arange(int(cx - w), int(cx + w))
Back   = int(cz - w), np.arange(int(cy - w), int(cy + w)), np.arange(int(cx - w), int(cx + w))

'''unit normal vectors at the perimeter of the bounding box'''
nxTop, nxBottom, nxRight, nxLeft, nxFront, nxBack = 0, 0, 1, -1, 0, 0
nyTop, nyBottom, nyRight, nyLeft, nyFront, nyBack = 1, -1, 0, 0, 0, 0
nzTop, nzBottom, nzRight, nzLeft, nzFront, nzBack = 0, 0, 0, 0, 1, -1


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




print(Fx)
print(Fy)
print(Fz)

print(Tx)
print(Ty)
print(Tz)




