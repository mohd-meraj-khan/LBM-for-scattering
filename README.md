# LBM for scattering of EM wave

## A parallel lattice Boltzmann solver for scattering and radiation force calculations

The solver employs the Lattice Boltzmann Method (LBM) as outlined by [Hauser and Verhey](https://doi.org/10.1103/PhysRevE.96.063306). 
To optimize performance, the LBM code is implemented in C and interfaced with Python using <b>ctypes</b>. Python was chosen for its ease of development and extensive library support. Additionally, OpenMP is utilized to parallelize the code for enhanced efficiency.

The solver is currently designed to compute radar cross-section (RCS) and radiation force for both 2D and 3D geometries under incident plane TM and TE waves. It operates effectively across all three scattering regimes: Rayleigh, Mie, and geometrical optics.







## Related publications

[Electromagnetic scattering by curved surfaces and calculation of radiation force: Lattice Boltzmann simulations](https://doi.org/10.1063/5.0234413)

*We compare LBM solutions with analytical solutions for smooth circular conducting and dielectric cylinders in scattering width and radiation force calculations. Additionally, we compare LBM solutions with [semi-analytical](https://doi.org/10.1364/OSAC.2.000298) solutions for corrugated elliptical conducting cylinders in radiation force calculations. In all cases, we find strong agreement between LBM solutions and both analytical and semi-analytical solutions across all three scattering regimes.*



