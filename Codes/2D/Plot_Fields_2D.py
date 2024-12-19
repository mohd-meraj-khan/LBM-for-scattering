import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
from Module_Parameters_2D import *


directory = 'data'
if not os.path.exists(directory):
    os.makedirs(directory)

plots = 'plots'
if not os.path.exists(plots):
    os.makedirs(plots)





'''loading the field data (cartesian)'''
Ex = np.load(directory+'/Ex_{}.npy'.format(ratio))
Ey = np.load(directory+'/Ey_{}.npy'.format(ratio))
Ez = np.load(directory+'/Ez_{}.npy'.format(ratio))

Hx = np.load(directory+'/Hx_{}.npy'.format(ratio))
Hy = np.load(directory+'/Hy_{}.npy'.format(ratio))
Hz = np.load(directory+'/Hz_{}.npy'.format(ratio))




###################################


plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)


################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(Ex), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/ExReal_a_{}.svg'.format(a))
plt.close()
##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(Ex), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/ExImag_a_{}.svg'.format(a))
plt.close()

##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(Ex), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/ExAbs_a_{}.svg'.format(a))
plt.close()

################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(Ey), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/EyReal_a_{}.svg'.format(a))
plt.close()
##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(Ey), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/EyImag_a_{}.svg'.format(a))
plt.close()
##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(Ey), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/EyAbs_a_{}.svg'.format(a))
plt.close()

################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(Ez), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/EzReal_a_{}.svg'.format(a))
plt.close()
##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(Ez), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/EzImag_a_{}.svg'.format(a))
plt.close()
##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(Ez), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/EzAbs_a_{}.svg'.format(a))
plt.close()
################################################################################


fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(Hx), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/HxReal_a_{}.svg'.format(a))
plt.close()
#####################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(Hx), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/HxImag_a_{}.svg'.format(a))
plt.close()
#####################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(Hx), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/HxAbs_a_{}.svg'.format(a))
plt.close()
################################################################################


fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(Hy), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/HyReal_a_{}.svg'.format(a))
plt.close()
######################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(Hy), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/HyImag_a_{}.svg'.format(a))
plt.close()
###################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(Hy), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/HyAbs_a_{}.svg'.format(a))
plt.close()

################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(Hz), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/HzReal_a_{}.svg'.format(a))
plt.close()
#####################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(Hz), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/HzImag_a_{}.svg'.format(a))
plt.close()
#####################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(Hz), cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.axhline(y=cy, linestyle='--', linewidth=0.5, color='k')
plt.axvline(x=cx, linestyle='--', linewidth=0.5, color='k')

plt.colorbar(shrink=0.97)

plt.xticks(np.linspace(0, Nx, 5))
plt.yticks(np.linspace(0, Ny, 5))


    
plt.savefig(plots+'/HzAbs_a_{}.svg'.format(a))
plt.close()
################################################################################

