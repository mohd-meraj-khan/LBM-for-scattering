import numpy as np


############################################# INSTANTANEOUS TRACTION FORCE ######################################

'''instantaneous traction vector in x direction'''
def fxInstantaneous(Ex, Ey, Ez, Hx, Hy, Hz, er1, mur1, nx, ny, nz):

    Txx = er1 * (Ex*Ex - 0.5 * (Ex*Ex + Ey*Ey + Ez*Ez)) + mur1 * (Hx*Hx - 0.5 * (Hx*Hx + Hy*Hy + Hz*Hz))
    Txy = er1 * (Ex*Ey) + mur1 * (Hx*Hy)
    Txz = er1 * (Ex*Ez) + mur1 * (Hx*Hz)

    fx = Txx * nx + Txy * ny + Txz * nz

    return fx



'''instantaneous traction vector in y direction'''
def fyInstantaneous(Ex, Ey, Ez, Hx, Hy, Hz, er1, mur1, nx, ny, nz):

    Tyx = er1 * (Ey*Ex) + mur1 * (Hy*Hx)
    Tyy = er1 * (Ey*Ey - 0.5 * (Ex*Ex + Ey*Ey + Ez*Ez)) + mur1 * (Hy*Hy - 0.5 * (Hx*Hx + Hy*Hy + Hz*Hz))
    Tyz = er1 * (Ey*Ez) + mur1 * (Hy*Hz)

    fy = Tyx * nx + Tyy * ny + Tyz * nz

    return fy



'''instantaneous traction vector in z direction'''
def fzInstantaneous(Ex, Ey, Ez, Hx, Hy, Hz, er1, mur1, nx, ny, nz):

    Tzx = er1 * (Ez*Ex) + mur1 * (Hz*Hx)
    Tzy = er1 * (Ez*Ey) + mur1 * (Hz*Hy)
    Tzz = er1 * (Ez*Ez - 0.5 * (Ex*Ex + Ey*Ey + Ez*Ez)) + mur1 * (Hz*Hz - 0.5 * (Hx*Hx + Hy*Hy + Hz*Hz))

    fz = Tzx * nx + Tzy * ny + Tzz * nz

    return fz

################################################################################################################




############################################# AVERAGE TRACTION FORCE ###########################################

'''average (over one time period) traction vector in x direction'''
def fxAverage(Ex, Ey, Ez, Hx, Hy, Hz, er1, mur1, nx, ny, nz):

    Ex_, Ey_, Ez_ = np.conjugate(Ex), np.conjugate(Ey), np.conjugate(Ez)
    Hx_, Hy_, Hz_ = np.conjugate(Hx), np.conjugate(Hy), np.conjugate(Hz)

    Txx = 0.5 * np.real(er1 * (Ex*Ex_ - 0.5 * (Ex*Ex_ + Ey*Ey_ + Ez*Ez_)) + mur1 * (Hx*Hx_ - 0.5 * (Hx*Hx_ + Hy*Hy_ + Hz*Hz_)))
    Txy = 0.5 * np.real(er1 * (Ex*Ey_) + mur1 * (Hx*Hy_))
    Txz = 0.5 * np.real(er1 * (Ex*Ez_) + mur1 * (Hx*Hz_))

    fx = Txx * nx + Txy * ny + Txz * nz

    return fx



'''average (over one time period) traction vector in y direction'''
def fyAverage(Ex, Ey, Ez, Hx, Hy, Hz, er1, mur1, nx, ny, nz):

    Ex_, Ey_, Ez_ = np.conjugate(Ex), np.conjugate(Ey), np.conjugate(Ez)
    Hx_, Hy_, Hz_ = np.conjugate(Hx), np.conjugate(Hy), np.conjugate(Hz)

    Tyx = 0.5 * np.real(er1 * (Ey*Ex_) + mur1 * (Hy*Hx_))
    Tyy = 0.5 * np.real(er1 * (Ey*Ey_ - 0.5 * (Ex*Ex_ + Ey*Ey_ + Ez*Ez_)) + mur1 * (Hy*Hy_ - 0.5 * (Hx*Hx_ + Hy*Hy_ + Hz*Hz_)))
    Tyz = 0.5 * np.real(er1 * (Ey*Ez_) + mur1 * (Hy*Hz_))

    fy = Tyx * nx + Tyy * ny + Tyz * nz

    return fy



'''average (over one time period) traction vector in z direction'''
def fzAverage(Ex, Ey, Ez, Hx, Hy, Hz, er1, mur1, nx, ny, nz):

    Ex_, Ey_, Ez_ = np.conjugate(Ex), np.conjugate(Ey), np.conjugate(Ez)
    Hx_, Hy_, Hz_ = np.conjugate(Hx), np.conjugate(Hy), np.conjugate(Hz)

    Tzx = 0.5 * np.real(er1 * (Ez*Ex_) + mur1 * (Hz*Hx_))
    Tzy = 0.5 * np.real(er1 * (Ez*Ey_) + mur1 * (Hz*Hy_))
    Tzz = 0.5 * np.real(er1 * (Ez*Ez_ - 0.5 * (Ex*Ex_ + Ey*Ey_ + Ez*Ez_)) + mur1 * (Hz*Hz_ - 0.5 * (Hx*Hx_ + Hy*Hy_ + Hz*Hz_)))

    fz = Tzx * nx + Tzy * ny + Tzz * nz

    return fz


#################################################################################################################





