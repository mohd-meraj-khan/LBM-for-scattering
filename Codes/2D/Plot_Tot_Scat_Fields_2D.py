import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
from Module_Parameters_2D import *


directory_scattered = 'data/scattered_field'
directory_total = 'data/total_field'
#####################

exreal = 'plots/plots_Field/real/ex'
if not os.path.exists(exreal):
    os.makedirs(exreal)

eyreal = 'plots/plots_Field/real/ey'
if not os.path.exists(eyreal):
    os.makedirs(eyreal)

ezreal = 'plots/plots_Field/real/ez'
if not os.path.exists(ezreal):
    os.makedirs(ezreal)

hxreal = 'plots/plots_Field/real/hx'
if not os.path.exists(hxreal):
    os.makedirs(hxreal)

hyreal = 'plots/plots_Field/real/hy'
if not os.path.exists(hyreal):
    os.makedirs(hyreal)

hzreal = 'plots/plots_Field/real/hz'
if not os.path.exists(hzreal):
    os.makedirs(hzreal)

########################

eximag = 'plots/plots_Field/imag/ex'
if not os.path.exists(eximag):
    os.makedirs(eximag)

eyimag = 'plots/plots_Field/imag/ey'
if not os.path.exists(eyimag):
    os.makedirs(eyimag)

ezimag = 'plots/plots_Field/imag/ez'
if not os.path.exists(ezimag):
    os.makedirs(ezimag)

hximag = 'plots/plots_Field/imag/hx'
if not os.path.exists(hximag):
    os.makedirs(hximag)

hyimag = 'plots/plots_Field/imag/hy'
if not os.path.exists(hyimag):
    os.makedirs(hyimag)

hzimag = 'plots/plots_Field/imag/hz'
if not os.path.exists(hzimag):
    os.makedirs(hzimag)

########################

exabs = 'plots/plots_Field/abs/ex'
if not os.path.exists(exabs):
    os.makedirs(exabs)

eyabs = 'plots/plots_Field/abs/ey'
if not os.path.exists(eyabs):
    os.makedirs(eyabs)

ezabs = 'plots/plots_Field/abs/ez'
if not os.path.exists(ezabs):
    os.makedirs(ezabs)

hxabs = 'plots/plots_Field/abs/hx'
if not os.path.exists(hxabs):
    os.makedirs(hxabs)

hyabs = 'plots/plots_Field/abs/hy'
if not os.path.exists(hyabs):
    os.makedirs(hyabs)

hzabs = 'plots/plots_Field/abs/hz'
if not os.path.exists(hzabs):
    os.makedirs(hzabs)

########################


'''loading the field data (cartesian)'''
ExScat = np.load(directory_scattered+'/ExScat_{}_{}.npy'.format(er2, ratio))
EyScat = np.load(directory_scattered+'/EyScat_{}_{}.npy'.format(er2, ratio))
EzScat = np.load(directory_scattered+'/EzScat_{}_{}.npy'.format(er2, ratio))
HxScat = np.load(directory_scattered+'/HxScat_{}_{}.npy'.format(er2, ratio))
HyScat = np.load(directory_scattered+'/HyScat_{}_{}.npy'.format(er2, ratio))
HzScat = np.load(directory_scattered+'/HzScat_{}_{}.npy'.format(er2, ratio))

ExScat[scatterer] = 0
EyScat[scatterer] = 0
EzScat[scatterer] = 0
HxScat[scatterer] = 0
HyScat[scatterer] = 0
HzScat[scatterer] = 0


ExTot = np.load(directory_total+'/ExTot_{}_{}.npy'.format(er2, ratio))
EyTot = np.load(directory_total+'/EyTot_{}_{}.npy'.format(er2, ratio))
EzTot = np.load(directory_total+'/EzTot_{}_{}.npy'.format(er2, ratio))
HxTot = np.load(directory_total+'/HxTot_{}_{}.npy'.format(er2, ratio))
HyTot = np.load(directory_total+'/HyTot_{}_{}.npy'.format(er2, ratio))
HzTot = np.load(directory_total+'/HzTot_{}_{}.npy'.format(er2, ratio))


###################################


def format_label(value, factor):
    normalized = value / factor
    # Check if the normalized value is an integer
    return f"{int(normalized)}" if normalized.is_integer() else f"{normalized:.1f}"




plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 2, lw = 0.75)
plt.rc('text', usetex = True)


################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(ExScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(exreal+'/ExScat_er_{}_Real_{}.svg'.format(er2, ratio))
plt.close()
##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(ExScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(eximag+'/ExScat_er_{}_Imag_{}.svg'.format(er2, ratio))
plt.close()

##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(ExScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(exabs+'/ExScat_er_{}_Abs_{}.svg'.format(er2, ratio))
plt.close()

################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(EyScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(eyreal+'/EyScat_er_{}_Real_{}.svg'.format(er2, ratio))
plt.close()
##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(EyScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(eyimag+'/EyScat_er_{}_Imag_{}.svg'.format(er2, ratio))
plt.close()
##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(EyScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(eyabs+'/EyScat_er_{}_Abs_{}.svg'.format(er2, ratio))
plt.close()

################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(EzScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(ezreal+'/EzScat_er_{}_Real_{}.svg'.format(er2, ratio))
plt.close()
##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(EzScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(ezimag+'/EzScat_er_{}_Imag_{}.svg'.format(er2, ratio))
plt.close()
##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(EzScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(ezabs+'/EzScat_er_{}_Abs_{}.svg'.format(er2, ratio))
plt.close()
################################################################################


fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(HxScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hxreal+'/HxScat_er_{}_Real_{}.svg'.format(er2, ratio))
plt.close()
#####################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(HxScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hximag+'/HxScat_er_{}_Imag_{}.svg'.format(er2, ratio))
plt.close()
#####################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(HxScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hxabs+'/HxScat_er_{}_Abs_{}.svg'.format(er2, ratio))
plt.close()
################################################################################


fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(HyScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hyreal+'/HyScat_er_{}_Real_{}.svg'.format(er2, ratio))
plt.close()
######################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(HyScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hyimag+'/HyScat_er_{}_Imag_{}.svg'.format(er2, ratio))
plt.close()
###################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(HyScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hyabs+'/HyScat_er_{}_Abs_{}.svg'.format(er2, ratio))
plt.close()

################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(HzScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hzreal+'/HzScat_er_{}_Real_{}.svg'.format(er2, ratio))
plt.close()
#####################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(HzScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hzimag+'/HzScat_er_{}_Imag_{}.svg'.format(er2, ratio))
plt.close()
#####################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(HzScat), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hzabs+'/HzScat_er_{}_Abs_{}.svg'.format(er2, ratio))
plt.close()
################################################################################




##########################################################################################################################################################

##########################################################################################################################################################



################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(ExTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(exreal+'/ExTot_er_{}_Real_{}.svg'.format(er2, ratio))
plt.close()
##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(ExTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(eximag+'/ExTot_er_{}_Imag_{}.svg'.format(er2, ratio))
plt.close()

##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(ExTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(exabs+'/ExTot_er_{}_Abs_{}.svg'.format(er2, ratio))
plt.close()

################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(EyTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(eyreal+'/EyTot_er_{}_Real_{}.svg'.format(er2, ratio))
plt.close()
##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(EyTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(eyimag+'/EyTot_er_{}_Imag_{}.svg'.format(er2, ratio))
plt.close()
##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(EyTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(eyabs+'/EyTot_er_{}_Abs_{}.svg'.format(er2, ratio))
plt.close()

################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(EzTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(ezreal+'/EzTot_er_{}_Real_{}.svg'.format(er2, ratio))
plt.close()
##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(EzTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(ezimag+'/EzTot_er_{}_Imag_{}.svg'.format(er2, ratio))
plt.close()
##################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(EzTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(ezabs+'/EzTot_er_{}_Abs_{}.svg'.format(er2, ratio))
plt.close()
################################################################################


fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(HxTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hxreal+'/HxTot_er_{}_Real_{}.svg'.format(er2, ratio))
plt.close()
#####################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(HxTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hximag+'/HxTot_er_{}_Imag_{}.svg'.format(er2, ratio))
plt.close()
#####################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(HxTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hxabs+'/HxTot_er_{}_Abs_{}.svg'.format(er2, ratio))
plt.close()
################################################################################


fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(HyTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hyreal+'/HyTot_er_{}_Real_{}.svg'.format(er2, ratio))
plt.close()
######################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(HyTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hyimag+'/HyTot_er_{}_Imag_{}.svg'.format(er2, ratio))
plt.close()
###################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(HyTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hyabs+'/HyTot_er_{}_Abs_{}.svg'.format(er2, ratio))
plt.close()

################################################################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.real(HzTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hzreal+'/HzTot_er_{}_Real_{}.svg'.format(er2, ratio))
plt.close()
#####################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.imag(HzTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hzimag+'/HzTot_er_{}_Imag_{}.svg'.format(er2, ratio))
plt.close()
#####################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize = (3.35, 2.75), dpi=600, constrained_layout = True)



plt.imshow(np.absolute(HzTot), vmin=-1, vmax=1, cmap='seismic', origin='lower')

circle1 = plt.Circle((cx, cy), a, linestyle='--', linewidth=0.5, color='k', fill=False)
plt.gca().add_patch(circle1)

plt.colorbar(shrink=0.97)

xticks = np.linspace(0, Nx, 5)
yticks = np.linspace(0, Ny, 5)

ax.set_xticks(xticks)
ax.set_xticklabels([format_label(x, a) for x in xticks])

ax.set_yticks(yticks)
ax.set_yticklabels([format_label(y, a) for y in xticks])

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')

    
plt.savefig(hzabs+'/HzTot_er_{}_Abs_{}.svg'.format(er2, ratio))
plt.close()
################################################################################



