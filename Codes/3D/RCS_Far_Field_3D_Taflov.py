import numpy as np
import scipy.special as sc
import sys
import os
from scipy.integrate import simpson
from Module_Parameters_3D import *



directory = 'data/scattered_field'


directory_rcs = 'data/rcs'
if not os.path.exists(directory_rcs):
    os.makedirs(directory_rcs)





'''loading the field data (cartesian)'''
Ex = np.load(directory+'/ExScat_{}_{}.npy'.format(er2, ratio))
Ey = np.load(directory+'/EyScat_{}_{}.npy'.format(er2, ratio))
Ez = np.load(directory+'/EzScat_{}_{}.npy'.format(er2, ratio))

Hx = np.load(directory+'/HxScat_{}_{}.npy'.format(er2, ratio))
Hy = np.load(directory+'/HyScat_{}_{}.npy'.format(er2, ratio))
Hz = np.load(directory+'/HzScat_{}_{}.npy'.format(er2, ratio))


'''unit normal vectors at the 6 faces of the bounding box'''
nTop = np.array([nxTop, nyTop, nzTop])
nBot = np.array([nxBot, nyBot, nzBot])
nRyt = np.array([nxRyt, nyRyt, nzRyt])
nLef = np.array([nxLef, nyLef, nzLef])
nFrt = np.array([nxFrt, nyFrt, nzFrt])
nBak = np.array([nxBak, nyBak, nzBak])


'''electric and magnetic fields at the 6 faces of the bounding box'''
E_Top = np.array([Ex[Top], Ey[Top], Ez[Top]])
E_Bot = np.array([Ex[Bot], Ey[Bot], Ez[Bot]])
E_Ryt = np.array([Ex[Ryt], Ey[Ryt], Ez[Ryt]])
E_Lef = np.array([Ex[Lef], Ey[Lef], Ez[Lef]])
E_Frt = np.array([Ex[Frt], Ey[Frt], Ez[Frt]])
E_Bak = np.array([Ex[Bak], Ey[Bak], Ez[Bak]])

H_Top = np.array([Hx[Top], Hy[Top], Hz[Top]])
H_Bot = np.array([Hx[Bot], Hy[Bot], Hz[Bot]])
H_Ryt = np.array([Hx[Ryt], Hy[Ryt], Hz[Ryt]])
H_Lef = np.array([Hx[Lef], Hy[Lef], Hz[Lef]])
H_Frt = np.array([Hx[Frt], Hy[Frt], Hz[Frt]])
H_Bak = np.array([Hx[Bak], Hy[Bak], Hz[Bak]])



'''electric and magnetic currents at the 6 faces of the bounding box'''
J_Top = np.cross(nTop, H_Top, axis=0)
J_Bot = np.cross(nBot, H_Bot, axis=0)
J_Ryt = np.cross(nRyt, H_Ryt, axis=0)
J_Lef = np.cross(nLef, H_Lef, axis=0)
J_Frt = np.cross(nFrt, H_Frt, axis=0)
J_Bak = np.cross(nBak, H_Bak, axis=0)

M_Top = - np.cross(nTop, E_Top, axis=0)
M_Bot = - np.cross(nBot, E_Bot, axis=0)
M_Ryt = - np.cross(nRyt, E_Ryt, axis=0)
M_Lef = - np.cross(nLef, E_Lef, axis=0)
M_Frt = - np.cross(nFrt, E_Frt, axis=0)
M_Bak = - np.cross(nBak, E_Bak, axis=0)


'''radii and angles of each grids of all the 6 faces of the bounding box'''
radTop, radBot, radRyt, radLef, radFrt, radBak = rad[Top], rad[Bot], rad[Ryt], rad[Lef], rad[Frt], rad[Bak]
theTop, theBot, theRyt, theLef, theFrt, theBak = the[Top], the[Bot], the[Ryt], the[Lef], the[Frt], the[Bak]
phiTop, phiBot, phiRyt, phiLef, phiFrt, phiBak = phi[Top], phi[Bot], phi[Ryt], phi[Lef], phi[Frt], phi[Bak]




'''grid for far-field'''
dthe = 1
dphi = 1

dthe_rad = np.deg2rad(dthe)
dphi_rad = np.deg2rad(dphi)



phii0 = phi0

the = np.arange(0, 360 + dthe, dthe)
phi = np.linspace(phii0 - dphi, phii0 + dphi, 3)

Phi_rad = np.deg2rad(phi)
The_rad = np.deg2rad(the)

Mesh_Phi, Mesh_The = np.meshgrid(Phi_rad, The_rad, indexing='ij')

The = np.ndarray.flatten(Mesh_The)
Phi = np.ndarray.flatten(Mesh_Phi)


'''initilization of vector potential at far distance'''
L_the, L_phi = [np.zeros_like(The, dtype=complex) for _ in range(2)]
N_the, N_phi = [np.zeros_like(The, dtype=complex) for _ in range(2)]



'''components of position vector at the fictitious boundary'''
x_dash_Top = radTop * np.sin(theTop) * np.cos(phiTop)
x_dash_Bot = radBot * np.sin(theBot) * np.cos(phiBot)
x_dash_Ryt = radRyt * np.sin(theRyt) * np.cos(phiRyt)
x_dash_Lef = radLef * np.sin(theLef) * np.cos(phiLef)
x_dash_Frt = radFrt * np.sin(theFrt) * np.cos(phiFrt)
x_dash_Bak = radBak * np.sin(theBak) * np.cos(phiBak)

y_dash_Top = radTop * np.sin(theTop) * np.sin(phiTop)
y_dash_Bot = radBot * np.sin(theBot) * np.sin(phiBot)
y_dash_Ryt = radRyt * np.sin(theRyt) * np.sin(phiRyt)
y_dash_Lef = radLef * np.sin(theLef) * np.sin(phiLef)
y_dash_Frt = radFrt * np.sin(theFrt) * np.sin(phiFrt)
y_dash_Bak = radBak * np.sin(theBak) * np.sin(phiBak)

z_dash_Top = radTop * np.cos(theTop)
z_dash_Bot = radBot * np.cos(theBot)
z_dash_Ryt = radRyt * np.cos(theRyt)
z_dash_Lef = radLef * np.cos(theLef)
z_dash_Frt = radFrt * np.cos(theFrt)
z_dash_Bak = radBak * np.cos(theBak)




'''calculation of vector potentials at far distance'''
for m in range(len(Phi)):

    ST = np.sin(The[m])
    SP = np.sin(Phi[m])
    CT = np.cos(The[m])
    CP = np.cos(Phi[m])

    r_top = x_dash_Top * ST * CP + y_dash_Top * ST * SP + z_dash_Top * CT
    r_bot = x_dash_Bot * ST * CP + y_dash_Bot * ST * SP + z_dash_Bot * CT
    r_Ryt = x_dash_Ryt * ST * CP + y_dash_Ryt * ST * SP + z_dash_Ryt * CT
    r_lef = x_dash_Lef * ST * CP + y_dash_Lef * ST * SP + z_dash_Lef * CT
    r_frt = x_dash_Frt * ST * CP + y_dash_Frt * ST * SP + z_dash_Frt * CT
    r_bak = x_dash_Bak * ST * CP + y_dash_Bak * ST * SP + z_dash_Bak * CT


    expTop, expBot = np.exp((0 + 1j)*k*r_top), np.exp((0 + 1j)*k*r_bot)
    expRyt, expLef = np.exp((0 + 1j)*k*r_Ryt), np.exp((0 + 1j)*k*r_lef)
    expFrt, expBak = np.exp((0 + 1j)*k*r_frt), np.exp((0 + 1j)*k*r_bak)



    L_the_top = simpson(simpson((M_Top[0] * CT * CP + M_Top[1] * CT * SP - M_Top[2] * ST) * expTop))
    L_the_bot = simpson(simpson((M_Bot[0] * CT * CP + M_Bot[1] * CT * SP - M_Bot[2] * ST) * expBot))
    L_the_Ryt = simpson(simpson((M_Ryt[0] * CT * CP + M_Ryt[1] * CT * SP - M_Ryt[2] * ST) * expRyt))
    L_the_lef = simpson(simpson((M_Lef[0] * CT * CP + M_Lef[1] * CT * SP - M_Lef[2] * ST) * expLef))
    L_the_frt = simpson(simpson((M_Frt[0] * CT * CP + M_Frt[1] * CT * SP - M_Frt[2] * ST) * expFrt))
    L_the_bak = simpson(simpson((M_Bak[0] * CT * CP + M_Bak[1] * CT * SP - M_Bak[2] * ST) * expBak))

    L_phi_top = simpson(simpson(( - M_Top[0] * SP + M_Top[1] * CP) * expTop))
    L_phi_bot = simpson(simpson(( - M_Bot[0] * SP + M_Bot[1] * CP) * expBot))
    L_phi_Ryt = simpson(simpson(( - M_Ryt[0] * SP + M_Ryt[1] * CP) * expRyt))
    L_phi_lef = simpson(simpson(( - M_Lef[0] * SP + M_Lef[1] * CP) * expLef))
    L_phi_frt = simpson(simpson(( - M_Frt[0] * SP + M_Frt[1] * CP) * expFrt))
    L_phi_bak = simpson(simpson(( - M_Bak[0] * SP + M_Bak[1] * CP) * expBak))




    N_the_top = simpson(simpson((J_Top[0] * CT * CP + J_Top[1] * CT * SP - J_Top[2] * ST) * expTop))
    N_the_bot = simpson(simpson((J_Bot[0] * CT * CP + J_Bot[1] * CT * SP - J_Bot[2] * ST) * expBot))
    N_the_Ryt = simpson(simpson((J_Ryt[0] * CT * CP + J_Ryt[1] * CT * SP - J_Ryt[2] * ST) * expRyt))
    N_the_lef = simpson(simpson((J_Lef[0] * CT * CP + J_Lef[1] * CT * SP - J_Lef[2] * ST) * expLef))
    N_the_frt = simpson(simpson((J_Frt[0] * CT * CP + J_Frt[1] * CT * SP - J_Frt[2] * ST) * expFrt))
    N_the_bak = simpson(simpson((J_Bak[0] * CT * CP + J_Bak[1] * CT * SP - J_Bak[2] * ST) * expBak))

    N_phi_top = simpson(simpson(( - J_Top[0] * SP + J_Top[1] * CP) * expTop))
    N_phi_bot = simpson(simpson(( - J_Bot[0] * SP + J_Bot[1] * CP) * expBot))
    N_phi_Ryt = simpson(simpson(( - J_Ryt[0] * SP + J_Ryt[1] * CP) * expRyt))
    N_phi_lef = simpson(simpson(( - J_Lef[0] * SP + J_Lef[1] * CP) * expLef))
    N_phi_frt = simpson(simpson(( - J_Frt[0] * SP + J_Frt[1] * CP) * expFrt))
    N_phi_bak = simpson(simpson(( - J_Bak[0] * SP + J_Bak[1] * CP) * expBak))

  
    
    

    L_the[m] = L_the_top + L_the_bot + L_the_Ryt + L_the_lef + L_the_frt + L_the_bak
    L_phi[m] = L_phi_top + L_phi_bot + L_phi_Ryt + L_phi_lef + L_phi_frt + L_phi_bak
    
    N_the[m] = N_the_top + N_the_bot + N_the_Ryt + N_the_lef + N_the_frt + N_the_bak
    N_phi[m] = N_phi_top + N_phi_bot + N_phi_Ryt + N_phi_lef + N_phi_frt + N_phi_bak


eta0 = 1    
    

L_the = L_the.reshape(Mesh_The.shape)
L_phi   = L_phi.reshape(Mesh_The.shape)

N_the = N_the.reshape(Mesh_The.shape)
N_phi   = N_phi.reshape(Mesh_The.shape)

'''far-field scattering width'''
RCS = k**2 / (4 * np.pi) * ( np.absolute(L_phi + eta0 * N_the)**2 + np.absolute(L_the - eta0 * N_phi)**2 ) / (np.pi * a**2)



np.save(directory_rcs+'/RCS_LBM_{}_{}_{}.npy'.format(ratio, er2, phii0), RCS[1])

##print(RCS[1])
