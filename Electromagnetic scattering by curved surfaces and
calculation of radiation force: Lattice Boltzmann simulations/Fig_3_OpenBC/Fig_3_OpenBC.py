import matplotlib.pyplot as plt
import numpy as np
import sys


nx = 6
ny = 6

cx = nx//2
cy = ny//2





plt.clf()

plt.rc('font', family = 'serif', size = 10)
plt.rc('xtick', labelsize = 10)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', markersize = 1.5, lw = 0.01)
plt.rc('text', usetex = True)

fig, ax = plt.subplots(figsize = (3.35, 3.35), dpi=600, constrained_layout = True)

# Plot the square grids
for i in range(nx):
    plt.plot([i, i], [0, nx], 'k--')

for j in range(ny):
    plt.plot([0, ny], [j, j], 'k--')






##################################

for i in range(nx+1):
    for j in range(ny+1):
        plt.scatter(i, j, color='k', marker='o')
##################################

# two headed arrow for delta x
plt.arrow(1, 4, 0, 1, head_width=0.06, head_length=0.15, linewidth=1, color='k', length_includes_head=True)
plt.arrow(1, 4, 0, -0.01, head_width=0.06, head_length=0.15, linewidth=1, color='k', length_includes_head=True)
plt.text(0.6, 4.5, r'$\Delta x$')

plt.arrow(1, 4, 1, 0, head_width=0.06, head_length=0.15, linewidth=1, color='k', length_includes_head=True)
plt.arrow(1, 4, -0.01, 0, head_width=0.06, head_length=0.15, linewidth=1, color='k', length_includes_head=True)
plt.text(1.35, 3.75, r'$\Delta x$')

# one head arrow
plt.arrow(nx, cy, -0.5, 0, head_width=0.06, head_length=0.15, linewidth=1, color='r', length_includes_head=True)
plt.arrow(nx-1, cy, -0.5, 0, head_width=0.06, head_length=0.15, linewidth=1, color='k', length_includes_head=True)

plt.text(nx-0.5, cy-0.5, r'$f_L^2$', color='r')
plt.text(nx-1-0.5, cy-0.5, r'$f_{L- \Delta x}^2$', color='k')


plt.arrow(cx, ny, 0, -0.5, head_width=0.06, head_length=0.15, linewidth=1, color='r', length_includes_head=True)
plt.arrow(cx, ny-1, 0, -0.5, head_width=0.06, head_length=0.15, linewidth=1, color='k', length_includes_head=True)

plt.text(cx+0.25, ny-0.5, r'$f_L^4$', color='r')
plt.text(cx+0.25, ny-1-0.5, r'$f_{L- \Delta x}^4$', color='k')

plt.arrow(cx, 0, 0, 0.5, head_width=0.06, head_length=0.15, linewidth=1, color='r', length_includes_head=True)
plt.arrow(cx, 1, 0, 0.5, head_width=0.06, head_length=0.15, linewidth=1, color='k', length_includes_head=True)


plt.arrow(1, 2, 0, 0.5, head_width=0.06, head_length=0.15, linewidth=1, color='g', length_includes_head=True)
plt.arrow(1, 2, 0, -0.5, head_width=0.06, head_length=0.15, linewidth=1, color='g', length_includes_head=True)
plt.arrow(1, 2, 0.5, 0, head_width=0.06, head_length=0.15, linewidth=1, color='g', length_includes_head=True)
plt.arrow(1, 2, -0.5, 0, head_width=0.06, head_length=0.15, linewidth=1, color='g', length_includes_head=True)

plt.text(1+0.05, 2+0.05, r'$f^0$', color='g')
plt.text(1+0.65, 2+0.0, r'$f^1$', color='g')
plt.text(0+0.25, 2+0.0, r'$f^2$', color='g')
plt.text(1+0.0, 2+0.65, r'$f^3$', color='g')
plt.text(1+0.0, 1+0.35, r'$f^4$', color='g')
plt.plot(1, 2, 'go')




##plt.text(22, 26, r'$a/ \Delta x = 6$', backgroundcolor='w')
##plt.text(22, 24, r'$L/a = 5$', backgroundcolor='w')



##plt.gca().add_patch(circle1)

plt.xlim(0, nx)
plt.ylim(0, ny)



plt.gca().set_aspect('equal', adjustable='box')
plt.gca().set_xlabel('')

plt.gca().tick_params(axis='both', direction='out')
plt.gca().set_xticks(np.linspace(0, nx, 7))
plt.gca().set_xticklabels([0, 1, 2, 3, 4, 5, 6])
plt.gca().set_yticks(np.linspace(0, ny, 7))
plt.gca().set_yticklabels([0, 1, 2, 3, 4, 5, 6])
plt.gca().set_xlabel(r'$x / \Delta x$')
plt.gca().set_ylabel(r'$y / \Delta x$')
##plt.gca().set_title(r'$a = {} \Delta x$'.format(a))



plt.savefig('openBC.svg')






