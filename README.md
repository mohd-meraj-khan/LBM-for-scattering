# LBM for scattering of EM wave

## A parallel lattice Boltzmann solver for scattering and radiation force calculations

The solver employs the Lattice Boltzmann Method (LBM) as outlined by [Hauser and Verhey](https://doi.org/10.1103/PhysRevE.96.063306). 
To optimize performance, the LBM code is implemented in C and interfaced with Python using <b>ctypes</b>. Python was chosen for its ease of development and extensive library support. Additionally, <b>OpenMP</b> is utilized to parallelize the code for enhanced efficiency.

The solver is currently designed to compute radar cross-section (RCS) and radiation force for both 2D and 3D geometries under incident plane TM and TE waves. It operates effectively across all three scattering regimes: Rayleigh, Mie, and geometrical optics.



> [!NOTE]
> I am currently seeking opportunities in the industry or a postdoc position in computational physics/biology. If you believe my skills and experience align with your needs, feel free to reach out to me at meraj@cacs.iitm.ac.in, connect with me on [LinkedIn](https://www.linkedin.com/in/meraj87), or view my [CV](https://drive.google.com/file/d/1gwV3vy8u4uV727nPQO8EAbnB6db8rSak/view?usp=sharing) for a quick overview. I look forward to connecting!


> [!NOTE]
> The code is currently a work in progress. The 2D version is stable and ready for use, while the 3D version requires further development to be fully functional.


> [!NOTE]
> If you have a research idea and are interested in collaborating, feel free to reach out to me at meraj@cacs.iitm.ac.in. I’d be happy to discuss potential opportunities!


## Related publications

[Electromagnetic scattering by curved surfaces and calculation of radiation force: Lattice Boltzmann simulations](https://doi.org/10.1063/5.0234413)

*We compare LBM solutions with analytical solutions for smooth circular conducting and dielectric cylinders in scattering width and radiation force calculations. Additionally, we compare LBM solutions with [semi-analytical](https://doi.org/10.1364/OSAC.2.000298) solutions for corrugated elliptical conducting cylinders in radiation force calculations. In all cases, we find strong agreement between LBM solutions and both analytical and semi-analytical solutions across all three scattering regimes.*

[Lattice Boltzmann method for electromagnetic wave scattering](https://arxiv.org/abs/2510.11042)

*In this paper, we propose the lattice Boltzmann method (LBM) as an alternative numerical approach for electromagnetic scattering. The method is systematically validated over a wide range of size parameters, thereby covering the Rayleigh, Mie, and geometric optics regimes, through comparison with established reference solutions. For circular cylinders, both perfect electrically conducting (PEC) and dielectric, LBM results are benchmarked against analytical Mie theory. For dielectric cylinders, comparisons are performed over a broad range of relative permittivities to assess accuracy across different material contrasts. Scattering from dielectric spheres is likewise compared with exact Mie solutions, showing excellent agreement. To assess performance for non-canonical geometries, we investigate a hexagonal dielectric cylinder and validate the results against the Discretized Mie-Formalism, demonstrating that LBM can accurately capture edge diffraction and sharp-facet effects. Overall, the study provides the first systematic benchmarking of LBM for electromagnetic scattering in one-, two-, and three-dimensional configurations, establishing it as a promising and versatile tool in computational electromagnetics.*

[Radiation Forces and Torques on Janus Cylinders](https://arxiv.org/abs/2509.22308)

*The interaction of electromagnetic waves with dielectric Janus particles gives rise to radiation forces and torques, governed by the dielectric properties, interface orientation, and the size-to-wavelength ratio. In this study, we employ the Lattice Boltzmann Method to compute the radiation-induced drag, lift, and torque on circular Janus cylinders when illuminated by a transverse magnetic polarized plane wave. We analyze both metallo-dielectric and dielectric Janus cylinders. For metallo-dielectric Janus cylinders, LBM predictions are validated against analytical results, showing excellent agreement in far-field bistatic scattering width, radiation force, and torque across a range of dielectric constants and interface orientations. Extending the study to dielectric Janus cylinders, we explore how the dielectric contrast and interface orientation shape the optomechanical response. Our findings show that radiation-induced forces and torques can be harnessed to drive and control the motion of dielectric Janus particles in optofluidic, active, and self-assembling systems. *




