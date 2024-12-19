import numpy as np
import scipy.special as sc
import sys
import os
from scipy.integrate import simpson
from Module_Parameters_2D import *



directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)




'''loading the field data (cartesian)'''
Ex = np.genfromtxt(directory+'/Ex_{}_{}.txt'.format(er2, theta), dtype=complex)
Ey = np.genfromtxt(directory+'/Ey_{}_{}.txt'.format(er2, theta), dtype=complex)
Ez = np.genfromtxt(directory+'/Ez_{}_{}.txt'.format(er2, theta), dtype=complex)

Hx = np.genfromtxt(directory+'/Hx_{}_{}.txt'.format(er2, theta), dtype=complex)
Hy = np.genfromtxt(directory+'/Hy_{}_{}.txt'.format(er2, theta), dtype=complex)
Hz = np.genfromtxt(directory+'/Hz_{}_{}.txt'.format(er2, theta), dtype=complex)



'''unit normal vectors at the 4 boundaries of the bounding box'''
nTop    = [nxTop,    nyTop,    nzTop]
nBottom = [nxBottom, nyBottom, nzBottom]
nRight  = [nxRight,  nyRight,  nzRight]
nLeft   = [nxLeft,   nyLeft,   nzLeft]
nFront  = [nxFront,  nyFront,  nzFront]
nBack   = [nxBack,   nyBack,   nzBack]



'''electric and magnetic fields at the 6 faces of the bounding box'''
E_Top    = [Ex[Top],    Ey[Top],    Ez[Top]]
E_Bottom = [Ex[Bottom], Ey[Bottom], Ez[Bottom]]
E_Right  = [Ex[Right],  Ey[Right],  Ez[Right]]
E_Left   = [Ex[Left],   Ey[Left],   Ez[Left]]
E_Front  = [Ex[Front],  Ey[Front],  Ez[Front]]
E_Back   = [Ex[Back],   Ey[Back],   Ez[Back]]

H_Top    = [Hx[Top],    Hy[Top],    Hz[Top]]
H_Bottom = [Hx[Bottom], Hy[Bottom], Hz[Bottom]]
H_Right  = [Hx[Right],  Hy[Right],  Hz[Right]]
H_Left   = [Hx[Left],   Hy[Left],   Hz[Left]]
H_Front  = [Hx[Front],  Hy[Front],  Hz[Front]]
H_Back   = [Hx[Back],   Hy[Back],   Hz[Back]]


'''electric and magnetic currents at the 4 boundaries of the bounding box'''
J_Top    =   np.cross(nTop,    H_Top,    axis=0)
J_Bottom =   np.cross(nBottom, H_Bottom, axis=0)
J_Right  =   np.cross(nRight,  H_Right,  axis=0)
J_Left   =   np.cross(nLeft,   H_Left,   axis=0)
J_Front  =   np.cross(nFront,  H_Front,  axis=0)
J_Back   =   np.cross(nBack,   H_Back,   axis=0)

M_Top    = - np.cross(nTop,    E_Top,    axis=0)
M_Bottom = - np.cross(nBottom, E_Bottom, axis=0)
M_Right  = - np.cross(nRight,  E_Right,  axis=0)
M_Left   = - np.cross(nLeft,   E_Left,   axis=0)
M_Front  = - np.cross(nFront,  E_Front,  axis=0)
M_Back   = - np.cross(nBack,   E_Back,   axis=0)


'''radii and angles of each grids of all the 4 boundaries of the bounding box'''
rTop, rBottom, rRight, rLeft, rFront, rBack = r[Top],   r[Bottom],   r[Right],   r[Left],   r[Front],   r[Back]
phiTop, phiBottom, phiRight, phiLeft, phiFront, phiBack = phi[Top], phi[Bottom], phi[Right], phi[Left], phi[Front], phi[Back]




'''grid for far-field'''
dr     = 1
dtheta = 1
dphi   = 1

dtheta_rad = np.deg2rad(dtheta)
dphi_rad   = np.deg2rad(dphi)

R0 = 100*a

phi0 = 0

r     = np.linspace(R0 - dr, R0 + dr, 3)
theta = np.arange(1, 179 + dtheta, dtheta)
phi   = np.linspace(phi0 - dphi, phi0 + dphi, 3)

Phi_rad   = np.deg2rad(phi)
Theta_rad = np.deg2rad(theta)

Mesh_Phi, Mesh_Theta, Mesh_R = np.meshgrid(Phi_rad, Theta_rad, r, indexing='ij')

R     = np.ndarray.flatten(Mesh_R)
Theta = np.ndarray.flatten(Mesh_Theta)
Phi   = np.ndarray.flatten(Mesh_Phi)


'''initilization of vector potential at far distance'''
A_r, A_theta, A_phi = [np.zeros_like(R, dtype=complex) for _ in range(3)]
F_r, F_theta, F_phi = [np.zeros_like(R, dtype=complex) for _ in range(3)]


'''components of position vector at the far distance'''
r_x = R * np.sin(Theta) * np.cos(Phi)
r_y = R * np.sin(Theta) * np.sin(Phi)
r_z = R * np.cos(Theta)


'''components of position vector at the fictitious boundary'''
r_dash_x_Top    = rTop    * np.sin(thetaTop)    * np.cos(phiTop)
r_dash_x_Bottom = rBottom * np.sin(thetaBottom) * np.cos(phiBottom)
r_dash_x_Right  = rRight  * np.sin(thetaRight)  * np.cos(phiRight)
r_dash_x_Left   = rLeft   * np.sin(thetaLeft)   * np.cos(phiLeft)
r_dash_x_Front  = rFront  * np.sin(thetaFront)  * np.cos(phiFront)
r_dash_x_Back   = rBack   * np.sin(thetaBack)   * np.cos(phiBack)

r_dash_y_Top    = rTop    * np.sin(thetaTop)    * np.sin(phiTop)
r_dash_y_Bottom = rBottom * np.sin(thetaBottom) * np.sin(phiBottom)
r_dash_y_Right  = rRight  * np.sin(thetaRight)  * np.sin(phiRight)
r_dash_y_Left   = rLeft   * np.sin(thetaLeft)   * np.sin(phiLeft)
r_dash_y_Front  = rFront  * np.sin(thetaFront)  * np.sin(phiFront)
r_dash_y_Back   = rBack   * np.sin(thetaBack)   * np.sin(phiBack)

r_dash_z_Top    = rTop    * np.cos(thetaTop)
r_dash_z_Bottom = rBottom * np.cos(thetaBottom)
r_dash_z_Right  = rRight  * np.cos(thetaRight)
r_dash_z_Left   = rLeft   * np.cos(thetaLeft)
r_dash_z_Front  = rFront  * np.cos(thetaFront)
r_dash_z_Back   = rBack   * np.cos(thetaBack)




'''calculation of vector potentials at far distance'''
for m in range(len(Phi)):

    distTop    = np.sqrt(np.square(r_x[m] - r_dash_x_Top)    + np.square(r_y[m] - r_dash_y_Top)    + np.square(r_z[m] - r_dash_z_Top))
    distBottom = np.sqrt(np.square(r_x[m] - r_dash_x_Bottom) + np.square(r_y[m] - r_dash_y_Bottom) + np.square(r_z[m] - r_dash_z_Bottom))
    distRight  = np.sqrt(np.square(r_x[m] - r_dash_x_Right)  + np.square(r_y[m] - r_dash_y_Right)  + np.square(r_z[m] - r_dash_z_Right))
    distLeft   = np.sqrt(np.square(r_x[m] - r_dash_x_Left)   + np.square(r_y[m] - r_dash_y_Left)   + np.square(r_z[m] - r_dash_z_Left))
    distFront  = np.sqrt(np.square(r_x[m] - r_dash_x_Front)  + np.square(r_y[m] - r_dash_y_Front)  + np.square(r_z[m] - r_dash_z_Front))
    distBack   = np.sqrt(np.square(r_x[m] - r_dash_x_Back)   + np.square(r_y[m] - r_dash_y_Back)   + np.square(r_z[m] - r_dash_z_Back))

    expTop,   expBottom = np.exp((0 - 1j)*k*distTop) / distTop,     np.exp((0 - 1j)*k*distBottom) / distBottom
    expRight, expLeft   = np.exp((0 - 1j)*k*distRight) / distRight, np.exp((0 - 1j)*k*distLeft) / distLeft
    expFront, expBack   = np.exp((0 - 1j)*k*distFront) / distFront, np.exp((0 - 1j)*k*distBack) / distBack


    Ax = simpson(J_Top[0]*expTop) + simpson(J_Bottom[0]*expBottom) + simpson(J_Right[0]*expRight) + simpson(J_Left[0]*expLeft + simpson(J_Front[0]*expFront) + simpson(J_Back[0]*expBack)
    Ay = simpson(J_Top[1]*expTop) + simpson(J_Bottom[1]*expBottom) + simpson(J_Right[1]*expRight) + simpson(J_Left[1]*expLeft + simpson(J_Front[1]*expFront) + simpson(J_Back[1]*expBack)
    Az = simpson(J_Top[2]*expTop) + simpson(J_Bottom[2]*expBottom) + simpson(J_Right[2]*expRight) + simpson(J_Left[2]*expLeft + simpson(J_Front[2]*expFront) + simpson(J_Back[2]*expBack)
    
    Fx = simpson(M_Top[0]*expTop) + simpson(M_Bottom[0]*expBottom) + simpson(M_Right[0]*expRight) + simpson(M_Left[0]*expLeft + simpson(M_Front[0]*expFront) + simpson(M_Back[0]*expBack)
    Fy = simpson(M_Top[1]*expTop) + simpson(M_Bottom[1]*expBottom) + simpson(M_Right[1]*expRight) + simpson(M_Left[1]*expLeft + simpson(M_Front[1]*expFront) + simpson(M_Back[1]*expBack)
    Fz = simpson(M_Top[2]*expTop) + simpson(M_Bottom[2]*expBottom) + simpson(M_Right[2]*expRight) + simpson(M_Left[2]*expLeft + simpson(M_Front[2]*expFront) + simpson(M_Back[2]*expBack)

    
    
    A_r[m]     = (0 - 1j) * 3 * mur1 / 4 * (   Ax * np.sin(Theta[m]) * np.cos(Phi[m]) + Ay * np.sin(Theta[m]) * np.sin(Phi[m]) + Az * np.cos(Theta[m]) )
    A_theta[m] = (0 - 1j) * 3 * mur1 / 4 * (   Ax * np.cos(Theta[m]) * np.cos(Phi[m]) + Ay * np.cos(Theta[m]) * np.sin(Phi[m]) - Az * np.sin(Theta[m]) )
    A_phi[m]   = (0 - 1j) * 3 * mur1 / 4 * ( - Ax * np.sin(Phi[m]) + Ay * np.cos(Phi[m]) )
    
    F_r[m]     = (0 - 1j) * 3 * er1  / 4 * (   Fx * np.sin(Theta[m]) * np.cos(Phi[m]) + Fy * np.sin(Theta[m]) * np.sin(Phi[m]) + Fz * np.cos(Theta[m]) )
    F_theta[m] = (0 - 1j) * 3 * er1  / 4 * (   Fx * np.cos(Theta[m]) * np.cos(Phi[m]) + Fy * np.cos(Theta[m]) * np.sin(Phi[m]) - Fz * np.sin(Theta[m]) )
    F_phi[m]   = (0 - 1j) * 3 * er1  / 4 * ( - Fx * np.sin(Phi[m]) + Fy * np.cos(Phi[m]) )
    
    

'''reshaping the vector potentials for gradient calculations'''
A_r     = A_r.reshape(Mesh_R.shape)
A_theta = A_theta.reshape(Mesh_R.shape)
A_phi   = A_phi.reshape(Mesh_R.shape)

F_r     = F_r.reshape(Mesh_R.shape)
F_theta = F_theta.reshape(Mesh_R.shape)
F_phi   = F_phi.reshape(Mesh_R.shape)


'''first order derivatives'''
d_Ar_dr     = np.gradient(A_r, dphi_rad, dtheta_rad, dr)[2]
d_Ar_dtheta = np.gradient(A_r, dphi_rad, dtheta_rad, dr)[1]
d_Ar_dphi   = np.gradient(A_r, dphi_rad, dtheta_rad, dr)[0]

d_Atheta_dr     = np.gradient(A_theta, dphi_rad, dtheta_rad, dr)[2]
d_Atheta_dtheta = np.gradient(A_theta, dphi_rad, dtheta_rad, dr)[1]
d_Atheta_dphi   = np.gradient(A_theta, dphi_rad, dtheta_rad, dr)[0]

d_Aphi_dr     = np.gradient(A_phi, dphi_rad, dtheta_rad, dr)[2]
d_Aphi_dtheta = np.gradient(A_phi, dphi_rad, dtheta_rad, dr)[1]
d_Aphi_dphi   = np.gradient(A_phi, dphi_rad, dtheta_rad, dr)[0]


d_Fr_dr     = np.gradient(F_r, dphi_rad, dtheta_rad, dr)[2]
d_Fr_dtheta = np.gradient(F_r, dphi_rad, dtheta_rad, dr)[1]
d_Fr_dphi   = np.gradient(F_r, dphi_rad, dtheta_rad, dr)[0]

d_Ftheta_dr     = np.gradient(F_theta, dphi_rad, dtheta_rad, dr)[2]
d_Ftheta_dtheta = np.gradient(F_theta, dphi_rad, dtheta_rad, dr)[1]
d_Ftheta_dphi   = np.gradient(F_theta, dphi_rad, dtheta_rad, dr)[0]

d_Fphi_dr     = np.gradient(F_phi, dphi_rad, dtheta_rad, dr)[2]
d_Fphi_dtheta = np.gradient(F_phi, dphi_rad, dtheta_rad, dr)[1]
d_Fphi_dphi   = np.gradient(F_phi, dphi_rad, dtheta_rad, dr)[0]




'''second order derivatives'''
d2_Ar_dr2         = np.gradient(d_Ar_dr, dphi_rad, dtheta_rad, dr)[2]
d2_Atheta_dtheta2 = np.gradient(d_Atheta_dtheta, dphi_rad, dtheta_rad, dr)[1]
d2_Aphi_dphi2     = np.gradient(d_Aphi_dphi, dphi_rad, dtheta_rad, dr)[0]

d2_Ar_dtheta_dr = np.gradient(d_Ar_dr, dphi_rad, dtheta_rad, dr)[1]
d2_Ar_dphi_dr   = np.gradient(d_Ar_dr, dphi_rad, dtheta_rad, dr)[0]

d2_Atheta_dr_dtheta   = np.gradient(d_Atheta_dtheta, dphi_rad, dtheta_rad, dr)[2]
d2_Atheta_dphi_dtheta = np.gradient(d_Atheta_dtheta, dphi_rad, dtheta_rad, dr)[0]

d2_Aphi_dr_dphi     = np.gradient(d_Aphi_dphi, dphi_rad, dtheta_rad, dr)[2]
d2_Aphi_dtheta_dphi = np.gradient(d_Aphi_dphi, dphi_rad, dtheta_rad, dr)[1]



'''gradient of divergence'''
TT = np.tan(Mesh_Theta)
ST = np.sin(Mesh_Theta)

grad_div_A_r     = d2_Ar_dr2 + 1 / Mesh_R * (2*d_Ar_dr + d_Atheta_dr/TT + d2_Atheta_dr_dtheta + d2_Aphi_dr_dphi/ST) - 1/Mesh_R**2*(2*A_r + A_theta/TT + d_Atheta_dtheta + d_Aphi_dphi/ST)
grad_div_A_theta = 1 / Mesh_R * d2_Ar_dtheta_dr + 1 / Mesh_R**2 * (2*d_Ar_dtheta - A_theta/ST**2 + d_Atheta_dtheta/TT + d2_Atheta_dtheta2 - d_Aphi_dphi/(ST*TT) + d2_Aphi_dtheta_dphi/ST)
grad_div_A_phi   = d2_Ar_dphi_dr/(Mesh_R*ST) + 1/Mesh_R**2 * (2*d_Ar_dphi/ST + d_Atheta_dphi/(ST*TT) + d2_Atheta_dphi_dtheta/ST + d2_Aphi_dphi2/ST**2)

grad_div_F_r     = d2_Fr_dr2 + 1 / Mesh_R * (2*d_Fr_dr + d_Ftheta_dr/TT + d2_Ftheta_dr_dtheta + d2_Fphi_dr_dphi/ST) - 1/Mesh_R**2*(2*F_r + F_theta/TT + d_Ftheta_dtheta + d_Fphi_dphi/ST)
grad_div_F_theta = 1 / Mesh_R * d2_Fr_dtheta_dr + 1 / Mesh_R**2 * (2*d_Fr_dtheta - F_theta/ST**2 + d_Ftheta_dtheta/TT + d2_Ftheta_dtheta2 - d_Fphi_dphi/(ST*TT) + d2_Fphi_dtheta_dphi/ST)
grad_div_F_phi   = d2_Fr_dphi_dr/(Mesh_R*ST) + 1/Mesh_R**2 * (2*d_Fr_dphi/ST + d_Ftheta_dphi/(ST*TT) + d2_Ftheta_dphi_dtheta/ST + d2_Fphi_dphi2/ST**2)


'''calculation of curl'''
curl_A_r     = A_phi/(Mesh_R*TT) + d_Aphi_dtheta/Mesh_R - d_Atheta_dphi/(Mesh_R*ST)
curl_A_theta = d_Ar_dphi/(Mesh_R*ST) - A_phi/Mesh_R - d_Aphi_dr
curl_A_phi   = A_theta/Mesh_R + d_Atheta_dr - d_Ar_dtheta/Mesh_R

curl_F_r     = F_phi/(Mesh_R*TT) + d_Fphi_dtheta/Mesh_R - d_Ftheta_dphi/(Mesh_R*ST)
curl_F_theta = d_Fr_dphi/(Mesh_R*ST) - F_phi/Mesh_R - d_Fphi_dr
curl_F_phi   = F_theta/Mesh_R + d_Ftheta_dr - d_Fr_dtheta/Mesh_R



'''far-field calculations'''
Er_far     = - (0 + 1j) * omega * (A_r     + 1 / k**2 * grad_div_A_r)     - 1 / (3 * er1)  * curl_F_r
Etheta_far = - (0 + 1j) * omega * (A_theta + 1 / k**2 * grad_div_A_theta) - 1 / (3 * er1)  * curl_F_theta
Ephi_far   = - (0 + 1j) * omega * (A_phi   + 1 / k**2 * grad_div_A_phi)   - 1 / (3 * er1)  * curl_F_phi

Hr_far     = - (0 + 1j) * omega * (F_r     + 1 / k**2 * grad_div_F_r)     + 1 / (3 * mur1) * curl_A_r
Htheta_far = - (0 + 1j) * omega * (F_theta + 1 / k**2 * grad_div_F_theta) + 1 / (3 * mur1) * curl_A_theta
Hphi_far   = - (0 + 1j) * omega * (F_phi   + 1 / k**2 * grad_div_F_phi)   + 1 / (3 * mur1) * curl_A_phi



Er_far     = np.absolute(Er_far)
Etheta_far = np.absolute(Etheta_far)
Ephi_far   = np.absolute(Ephi_far)

Hr_far     = np.absolute(Hr_far)
Htheta_far = np.absolute(Htheta_far)
Hphi_far   = np.absolute(Hphi_far)


'''far-field scattering width'''
RCS = 4 * (R0/a)**2 * (Er_far[1, :, 1]**2 + Etheta_far[1, :, 1]**2 + Ephi_far[1, :, 1]**2)



brcsFar = open(directory+"/RCS_LBM_{}.txt".format(er2), "w")
np.savetxt(brcsFar, RCS)
brcsFar.close()


