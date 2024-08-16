


############################################# TRACTION FORCE ###################################################

############################### traction vector fx function ###################
def tractionFx(Ex, Ey, Ez, Hx, Hy, Hz, er1, mur1, nx, ny, nz):

    fx = 0

    Txx = 0
    Txx = er1 * (Ex*Ex - 0.5 * (Ex*Ex + Ey*Ey + Ez*Ez)) + mur1 * (Hx*Hx - 0.5 * (Hx*Hx + Hy*Hy + Hz*Hz))

    Txy = 0
    Txy = er1 * (Ex*Ey) + mur1 * (Hx*Hy)

    Txz = 0
    Txz = er1 * (Ex*Ez) + mur1 * (Hx*Hz)

    fx = Txx * nx + Txy * ny + Txz * nz

    return fx
###############################################################################


############################### traction vector fy function ###################
def tractionFy(Ex, Ey, Ez, Hx, Hy, Hz, er1, mur1, nx, ny, nz):

    fy = 0

    Tyx = 0
    Tyx = er1 * (Ey*Ex) + mur1 * (Hy*Hx)

    Tyy = 0
    Tyy = er1 * (Ey*Ey - 0.5 * (Ex*Ex + Ey*Ey + Ez*Ez)) + mur1 * (Hy*Hy - 0.5 * (Hx*Hx + Hy*Hy + Hz*Hz))

    Tyz = 0
    Tyz = er1 * (Ey*Ez) + mur1 * (Hy*Hz)

    fy = Tyx * nx + Tyy * ny + Tyz * nz

    return fy
###############################################################################


############################### traction vector fz function ###################
def tractionFz(Ex, Ey, Ez, Hx, Hy, Hz, er1, mur1, nx, ny, nz):

    fz = 0

    Tzx = 0
    Tzx = er1 * (Ez*Ex) + mur1 * (Hz*Hx)

    Tzy = 0
    Tzy = er1 * (Ez*Ey) + mur1 * (Hz*Hy)

    Tzz = 0
    Tzz = er1 * (Ez*Ez - 0.5 * (Ex*Ex + Ey*Ey + Ez*Ez)) + mur1 * (Hz*Hz - 0.5 * (Hx*Hx + Hy*Hy + Hz*Hz))

    fz = Tzx * nx + Tzy * ny + Tzz * nz

    return fz
###############################################################################

###########################################################################################################################





