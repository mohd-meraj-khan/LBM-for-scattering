import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import matplotlib as mpl
from Module_Parameters_3D import *


pv.global_theme.font.family = "times"
pv.global_theme.font.size = 25

# --------------------------
# Load field data
# --------------------------
Ex = np.load(f'data/total_field/ExTot_{er2}_{ratio}.npy')
Ey = np.load(f'data/total_field/EyTot_{er2}_{ratio}.npy')
Ez = np.load(f'data/total_field/EzTot_{er2}_{ratio}.npy')

# Scalar field (magnitude of E-field)
scalar_field = np.sqrt(np.abs(Ex)**2 + np.abs(Ey)**2 + np.abs(Ez)**2)

nz, ny, nx = scalar_field.shape

# Coordinates scaled by a
x = np.arange(nx) / a
y = np.arange(ny) / a
z = np.arange(nz) / a
Z, Y, X = np.meshgrid(z, y, x, indexing='ij')

# Create structured grid
grid = pv.StructuredGrid()
grid.points = np.column_stack((Z.ravel(), Y.ravel(), X.ravel()))
grid.dimensions = (nz, ny, nx)
grid["E_mag"] = scalar_field.flatten(order="F")

# Sphere (scatterer)
cx, cy, cz = (nx - 1) / (2*a), (ny - 1) / (2*a), (nz - 1) / (2*a)
sphere = pv.Sphere(radius=1, center=(cx, cy, cz), theta_resolution=64, phi_resolution=64)

# --------------------------
# PyVista render (NO scalar bar)
# --------------------------

dpi = 600
width_in, height_in = 3.35, 2.5
window_size = (int(width_in * dpi), int(height_in * dpi))

values = np.linspace(np.min(scalar_field), np.max(scalar_field), 5)

plotter = pv.Plotter(off_screen=True, window_size=window_size)
contours = grid.contour(isosurfaces=values, scalars="E_mag")

plotter.add_mesh(contours, cmap="cool", opacity=0.4, show_scalar_bar=False)  # 🔴 no scalar bar
plotter.add_mesh(sphere, color="grey", opacity=0.4, specular=1, specular_power=100, smooth_shading=True)

##plotter.show_grid(color="black", xtitle="x/a", ytitle="y/a", ztitle="z/a")

##plotter.add_axes()


##plotter.view_vector((-1, 0, 0), (1, 0, 0))


# Zooming the 3D render to fit in the white space
zoom_factor = 1.3  # Increase this value (e.g., 2.0) to zoom in more, decrease (e.g., 0.5) to zoom out
current_position = plotter.camera.position
focal_point = plotter.camera.focal_point
distance = np.linalg.norm(np.array(current_position) - np.array(focal_point))
new_distance = distance / zoom_factor
direction = np.array(current_position) - np.array(focal_point)
new_position = np.array(focal_point) + (direction / np.linalg.norm(direction)) * new_distance
plotter.camera.position = tuple(new_position)

# Screenshot (just the 3D render)
img = plotter.screenshot(return_img=True)



# --------------------------
# Matplotlib plot
# --------------------------




plt.rc('font', family='serif', size=10)
plt.rc('xtick', labelsize=10)
plt.rc('ytick', labelsize=10)
plt.rc('lines', markersize=2, lw=0.75)
plt.rc('text', usetex=True)

fig, ax = plt.subplots(figsize=(width_in, height_in), dpi=dpi)
ax.imshow(img, extent=[0, width_in*dpi, 0, height_in*dpi])  # Match image to figure size


ax.axis("off")


### Draw X-axis arrow
##ax.annotate("", xy=(Nx*1, 0), xytext=(0, 0),
##            arrowprops=dict(arrowstyle='->', color="k", lw=0.5))
##
### Draw Y-axis arrow
##ax.annotate("", xy=(0, Ny*1), xytext=(0, 0),
##            arrowprops=dict(arrowstyle='->', color="k", lw=0.5))
##
### Axis labels
##ax.text(1.1*Nx, 0.1*Ny, "z", ha="center", va="top")
##ax.text(0.1*Nx, 1.1*Ny, "y", ha="right", va="center")




# Create Matplotlib colorbar
norm = mpl.colors.Normalize(vmin=scalar_field.min(), vmax=scalar_field.max())
cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="cool"), ax=ax, shrink=0.6, aspect=30)
cbar.set_label(r"$|{\bf E}^{tot}|/E_0$")







plt.savefig("sphere_er_{}_ratio_{}.svg".format(er2, ratio), dpi=dpi, bbox_inches="tight")
