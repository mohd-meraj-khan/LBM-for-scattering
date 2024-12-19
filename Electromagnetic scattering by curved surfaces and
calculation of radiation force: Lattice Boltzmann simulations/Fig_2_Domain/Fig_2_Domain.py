import matplotlib.pyplot as plt
import numpy as np
import sys




# domain size
nx = 30
ny = 30

# center of scatterer
cx = nx//2 + 0.5
cy = ny//2 + 0.5

# size of scatterer
a = 6


ratio = 2
wavelength = a / ratio
k = 2 * np.pi / wavelength


# sinusoidal wave
x = np.arange(0, nx, 0.001)
y = np.zeros(len(x))
y = 4*np.sin(k*x) + ny//2 + 0.5


# cartesian to polar conversion
r = np.zeros([ny+1, nx+1])
phi = np.zeros([ny+1, nx+1])

for i in range(ny+1):
    for j in range(nx+1):
        r[i,j] = np.sqrt((i-ny//2 - 0.5)**2 + (j-nx//2 - 0.5)**2)
        phi[i,j] = np.arctan2((i-ny//2 - 0.5),(j-nx//2 - 0.5))





plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 0.75, lw = 0.01)
plt.rc('text', usetex = True)


fig, ax = plt.subplots(figsize = (3.35, 3.35), dpi=600, constrained_layout = True)


# Plot the square grids
for i in range(nx):
    plt.plot([i, i], [0, nx], 'k--')

for j in range(ny):
    plt.plot([0, ny], [j, j], 'k--')



circle1 = plt.Circle((cx, cy), a, color='k', fill=False)


# ploting nodes
for i in range(nx+1):
    for j in range(ny+1):
        if r[i, j] <= a :
            plt.scatter(i, j, color='r', marker='o')
        else:
            plt.scatter(i, j, color='k', marker='o')



# wave plot
plt.plot(x[:6000], y[:6000], 'r-', lw=0.75)

# centerlines plot
plt.axhline(y=cy, color='k', lw=0.5, linestyle='--')
plt.axvline(x=cx, color='k', lw=0.5, linestyle='--')


# vertical lines for wavelength measurment
plt.axvline(x=3, ymin=0.5, ymax=0.7, lw=0.5, linestyle='--', color='k')
plt.axvline(x=6, ymin=0.5, ymax=0.7, lw=0.5, linestyle='--', color='k')


# two headed arrow for lambda
plt.arrow(3, 20.5, 3, 0, head_width=0.2, head_length=0.5, linewidth=0.5, color='k', length_includes_head=True)
plt.arrow(3.1, 20.5, -0.1, 0, head_width=0.2, head_length=0.5, linewidth=0.5, color='k', length_includes_head=True)
plt.text(4.21, 20.65, r'$\lambda$')


# one head arrow for radius
plt.arrow(cx, cy, a*np.cos(0.5), a*np.sin(0.5), head_width=0.3, head_length=0.5, linewidth=0.5, color='k', length_includes_head=True)
plt.text(cx+2, cy+2, r'$a$')




plt.gca().add_patch(circle1)

plt.xlim(0, nx)
plt.ylim(0, ny)



plt.gca().set_aspect('equal', adjustable='box')
plt.gca().set_xlabel('')



plt.gca().tick_params(axis='both', direction='out')
plt.gca().set_xticks(np.linspace(0, nx, 6))
plt.gca().set_xticklabels([0, 1, 2, 3, 4, 5])
plt.gca().set_yticks(np.linspace(0, ny, 6))
plt.gca().set_yticklabels([0, 1, 2, 3, 4, 5])
plt.gca().set_xlabel(r'$x/a$')
plt.gca().set_ylabel(r'$y/a$')
##plt.gca().set_title(r'$a = {} \Delta x$'.format(a))



plt.savefig('DomainSchematicA_1.0Smooth.svg')






