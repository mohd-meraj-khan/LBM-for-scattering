import numpy as np
import scipy.special as sc
import sys
import os
from scipy.integrate import simpson
from Module_Parameters_2D import *



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




'''unit normal vectors at the 4 boundaries of the bounding box'''
nTop    = [nxTop,    nyTop,    nzTop]
nRight  = [nxRight,  nyRight,  nzRight]
nBottom = [nxBottom, nyBottom, nzBottom]
nLeft   = [nxLeft,   nyLeft,   nzLeft]



'''electric and magnetic fields at the 4 boundaries of the bounding box'''
E_Top    = [Ex[Top],    Ey[Top],    Ez[Top]]
E_Right  = [Ex[Right],  Ey[Right],  Ez[Right]]
E_Bottom = [Ex[Bottom], Ey[Bottom], Ez[Bottom]]
E_Left   = [Ex[Left],   Ey[Left],   Ez[Left]]

H_Top    = [Hx[Top],    Hy[Top],    Hz[Top]]
H_Right  = [Hx[Right],  Hy[Right],  Hz[Right]]
H_Bottom = [Hx[Bottom], Hy[Bottom], Hz[Bottom]]
H_Left   = [Hx[Left],   Hy[Left],   Hz[Left]]


'''electric and magnetic currents at the 4 boundaries of the bounding box'''
J_Top    =   np.cross(nTop,    H_Top,    axis=0)
J_Right  =   np.cross(nRight,  H_Right,  axis=0)
J_Bottom =   np.cross(nBottom, H_Bottom, axis=0)
J_Left   =   np.cross(nLeft,   H_Left,   axis=0)

M_Top    = - np.cross(nTop,    E_Top,    axis=0)
M_Right  = - np.cross(nRight,  E_Right,  axis=0)
M_Bottom = - np.cross(nBottom, E_Bottom, axis=0)
M_Left   = - np.cross(nLeft,   E_Left,   axis=0)


'''radii and angles of each grids of all the 4 boundaries of the bounding box'''
rhoTop, rhoRight, rhoBottom, rhoLeft = r[Top],   r[Right],   r[Bottom],   r[Left]
phiTop, phiRight, phiBottom, phiLeft = phi[Top], phi[Right], phi[Bottom], phi[Left]




'''grid for far-field'''
drho = 1
dphi = 1

dphi_rad = np.deg2rad(dphi)

rho0 = 100*a

rho = np.linspace(rho0 - drho, rho0 + drho, 3)
phi = np.arange(0, 360 + dphi, dphi)

Phi_rad = np.deg2rad(phi)

Mesh_Phi, Mesh_Rho = np.meshgrid(Phi_rad, rho, indexing='ij')

Rho = np.ndarray.flatten(Mesh_Rho)
Phi = np.ndarray.flatten(Mesh_Phi)


'''initilization of vector potential at far distance'''
A_r, A_theta, A_z = [np.zeros_like(Rho, dtype=complex) for _ in range(3)]
F_r, F_theta, F_z = [np.zeros_like(Rho, dtype=complex) for _ in range(3)]


'''components of position vector at the far distance'''
rho_x = Rho * np.cos(Phi)
rho_y = Rho * np.sin(Phi)
rho_z = np.zeros_like(rho_x)


'''components of position vector at the fictitious boundary'''
rho_dash_x_Top    = rhoTop    * np.cos(phiTop)
rho_dash_x_Right  = rhoRight  * np.cos(phiRight)
rho_dash_x_Bottom = rhoBottom * np.cos(phiBottom)
rho_dash_x_Left   = rhoLeft   * np.cos(phiLeft)

rho_dash_y_Top    = rhoTop    * np.sin(phiTop)
rho_dash_y_Right  = rhoRight  * np.sin(phiRight)
rho_dash_y_Bottom = rhoBottom * np.sin(phiBottom)
rho_dash_y_Left   = rhoLeft   * np.sin(phiLeft)

rho_dash_z_Top    = np.zeros_like(rho_dash_x_Top)
rho_dash_z_Right  = np.zeros_like(rho_dash_x_Right)
rho_dash_z_Bottom = np.zeros_like(rho_dash_x_Bottom)
rho_dash_z_Left   = np.zeros_like(rho_dash_x_Left)




'''calculation of vector potentials at far distance'''
for m in range(len(Phi)):

    distTop    = np.sqrt(np.square(rho_x[m] - rho_dash_x_Top)    + np.square(rho_y[m] - rho_dash_y_Top)    + np.square(rho_z[m] - rho_dash_z_Top))
    distRight  = np.sqrt(np.square(rho_x[m] - rho_dash_x_Right)  + np.square(rho_y[m] - rho_dash_y_Right)  + np.square(rho_z[m] - rho_dash_z_Right))
    distBottom = np.sqrt(np.square(rho_x[m] - rho_dash_x_Bottom) + np.square(rho_y[m] - rho_dash_y_Bottom) + np.square(rho_z[m] - rho_dash_z_Bottom))
    distLeft   = np.sqrt(np.square(rho_x[m] - rho_dash_x_Left)   + np.square(rho_y[m] - rho_dash_y_Left)   + np.square(rho_z[m] - rho_dash_z_Left))

    HankelTop, HankelRight, HankelBottom, HankelLeft = sc.hankel2(0, k*distTop), sc.hankel2(0, k*distRight), sc.hankel2(0, k*distBottom), sc.hankel2(0, k*distLeft)


    Ax = simpson(J_Top[0]*HankelTop) + simpson(J_Right[0]*HankelRight) + simpson(J_Bottom[0]*HankelBottom) + simpson(J_Left[0]*HankelLeft)
    Ay = simpson(J_Top[1]*HankelTop) + simpson(J_Right[1]*HankelRight) + simpson(J_Bottom[1]*HankelBottom) + simpson(J_Left[1]*HankelLeft)
    Az = simpson(J_Top[2]*HankelTop) + simpson(J_Right[2]*HankelRight) + simpson(J_Bottom[2]*HankelBottom) + simpson(J_Left[2]*HankelLeft)
    
    Fx = simpson(M_Top[0]*HankelTop) + simpson(M_Right[0]*HankelRight) + simpson(M_Bottom[0]*HankelBottom) + simpson(M_Left[0]*HankelLeft)
    Fy = simpson(M_Top[1]*HankelTop) + simpson(M_Right[1]*HankelRight) + simpson(M_Bottom[1]*HankelBottom) + simpson(M_Left[1]*HankelLeft)
    Fz = simpson(M_Top[2]*HankelTop) + simpson(M_Right[2]*HankelRight) + simpson(M_Bottom[2]*HankelBottom) + simpson(M_Left[2]*HankelLeft)

    
    
    A_r[m]     = (0 - 1j) * 3 * mur1 / 4 * (   Ax * np.cos(Phi[m]) + Ay * np.sin(Phi[m]) )
    A_theta[m] = (0 - 1j) * 3 * mur1 / 4 * ( - Ax * np.sin(Phi[m]) + Ay * np.cos(Phi[m]) )
    A_z[m]     = (0 - 1j) * 3 * mur1 / 4 * Az
    
    F_r[m]     = (0 - 1j) * 3 * er1  / 4 * (   Fx * np.cos(Phi[m]) + Fy * np.sin(Phi[m]) )
    F_theta[m] = (0 - 1j) * 3 * er1  / 4 * ( - Fx * np.sin(Phi[m]) + Fy * np.cos(Phi[m]) )
    F_z[m]     = (0 - 1j) * 3 * er1  / 4 * Fz
    
    

'''reshaping the vector potentials for gradient calculations'''
A_r     = A_r.reshape(Mesh_Rho.shape)
A_theta = A_theta.reshape(Mesh_Rho.shape)
A_z     = A_z.reshape(Mesh_Rho.shape)

F_r     = F_r.reshape(Mesh_Rho.shape)
F_theta = F_theta.reshape(Mesh_Rho.shape)
F_z     = F_z.reshape(Mesh_Rho.shape)


'''first order derivatives'''
del_Ar_del_r     = np.gradient(A_r, dphi_rad, drho)[1]
del_Ar_del_phi   = np.gradient(A_r, dphi_rad, drho)[0]

del_Aphi_del_r   = np.gradient(A_theta, dphi_rad, drho)[1]
del_Aphi_del_phi = np.gradient(A_theta, dphi_rad, drho)[0]

del_Az_del_r     = np.gradient(A_z, dphi_rad, drho)[1]
del_Az_del_phi   = np.gradient(A_z, dphi_rad, drho)[0]


del_Fr_del_r     = np.gradient(F_r, dphi_rad, drho)[1]
del_Fr_del_phi   = np.gradient(F_r, dphi_rad, drho)[0]

del_Fphi_del_r   = np.gradient(F_theta, dphi_rad, drho)[1]
del_Fphi_del_phi = np.gradient(F_theta, dphi_rad, drho)[0]

del_Fz_del_r     = np.gradient(F_z, dphi_rad, drho)[1]
del_Fz_del_phi   = np.gradient(F_z, dphi_rad, drho)[0]



'''second order derivatives'''
del2_Ar_del_r2          = np.gradient(del_Ar_del_r,     dphi_rad, drho)[1]
del2_Aphi_del_phi2      = np.gradient(del_Aphi_del_phi, dphi_rad, drho)[0]
del2_Ar_del_phi_del_r   = np.gradient(del_Ar_del_r,     dphi_rad, drho)[0]
del2_Aphi_del_r_del_phi = np.gradient(del_Aphi_del_phi, dphi_rad, drho)[1]

del2_Fr_del_r2          = np.gradient(del_Fr_del_r,     dphi_rad, drho)[1]
del2_Fphi_del_phi2      = np.gradient(del_Fphi_del_phi, dphi_rad, drho)[0]
del2_Fr_del_phi_del_r   = np.gradient(del_Fr_del_r,     dphi_rad, drho)[0]
del2_Fphi_del_r_del_phi = np.gradient(del_Fphi_del_phi, dphi_rad, drho)[1]


'''gradient of divergence'''
grad_div_A_r   = del2_Ar_del_r2 + 1 / Mesh_Rho * (del_Ar_del_r + del2_Aphi_del_r_del_phi) - 1 / Mesh_Rho**2 * (A_r + del_Aphi_del_phi) 
grad_div_A_phi = 1 / Mesh_Rho * del2_Ar_del_phi_del_r + 1 / Mesh_Rho**2 * (del_Ar_del_phi +  del2_Aphi_del_phi2)
grad_div_A_z   = np.zeros_like(grad_div_A_r)

grad_div_F_r   = del2_Fr_del_r2 + 1 / Mesh_Rho * (del_Fr_del_r + del2_Fphi_del_r_del_phi) - 1 / Mesh_Rho**2 * (F_r + del_Fphi_del_phi) 
grad_div_F_phi = 1 / Mesh_Rho * del2_Fr_del_phi_del_r + 1 / Mesh_Rho**2 * (del_Fr_del_phi +  del2_Fphi_del_phi2)
grad_div_F_z   = np.zeros_like(grad_div_F_r)

'''calculation of curl'''
curl_A_r   =   del_Az_del_phi / Mesh_Rho
curl_A_phi = - del_Az_del_r
curl_A_z   =   del_Aphi_del_r - del_Ar_del_phi / Mesh_Rho

curl_F_r   =   del_Fz_del_phi / Mesh_Rho
curl_F_phi = - del_Fz_del_r
curl_F_z   =   del_Fphi_del_r - del_Fr_del_phi / Mesh_Rho


'''far-field calculations'''
Er_far   = - (0 + 1j) * omega * (A_r     + 1 / k**2 * grad_div_A_r)   - 1 / (3 * er1) * curl_F_r
Ephi_far = - (0 + 1j) * omega * (A_theta + 1 / k**2 * grad_div_A_phi) - 1 / (3 * er1) * curl_F_phi
Ez_far   = - (0 + 1j) * omega * (A_z     + 1 / k**2 * grad_div_A_z)   - 1 / (3 * er1) * curl_F_z

Hr_far   = - (0 + 1j) * omega * (F_r     + 1 / k**2 * grad_div_F_r)   + 1 / (3 * mur1) * curl_A_r
Hphi_far = - (0 + 1j) * omega * (F_theta + 1 / k**2 * grad_div_F_phi) + 1 / (3 * mur1) * curl_A_phi
Hz_far   = - (0 + 1j) * omega * (F_z     + 1 / k**2 * grad_div_F_z)   + 1 / (3 * mur1) * curl_A_z


Er_far   = np.absolute(Er_far)
Ephi_far = np.absolute(Ephi_far)
Ez_far   = np.absolute(Ez_far)

Hr_far   = np.absolute(Hr_far)
Hphi_far = np.absolute(Hphi_far)
Hz_far   = np.absolute(Hz_far)


'''far-field scattering width'''
RCS = 2*np.pi*rho0/wavelength*Ez_far[:, 1]**2



np.save(directory_rcs+'/RCS_LBM_{}_{}.npy'.format(er2, ratio), RCS)



