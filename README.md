# LBM for scattering of EM wave

## A parallel lattice Boltzmann solver for scattering and radiation force calculations

The solver employs the Lattice Boltzmann Method (LBM) as outlined by [Hauser and Verhey](https://doi.org/10.1103/PhysRevE.96.063306). 
To optimize performance, the LBM code is implemented in C and interfaced with Python using <b>ctypes</b>. Python was chosen for its ease of development and extensive library support. Additionally, <b>OpenMP</b> is utilized to parallelize the code for enhanced efficiency.

The solver is currently designed to compute radar cross-section (RCS) and radiation force for both 2D and 3D geometries under incident plane TM and TE waves. It operates effectively across all three scattering regimes: Rayleigh, Mie, and geometrical optics.



> [!NOTE]
> If you have a research idea and are interested in collaborating, feel free to reach out to me at meraj@cacs.iitm.ac.in. I’d be happy to discuss potential opportunities!

## 📘 Thesis


### [**A Lattice Boltzmann Framework for Simulating Electromagnetic Scattering and Radiation Forces**](https://doi.org/10.5281/zenodo.19487120)

A unified Lattice Boltzmann framework for electromagnetic scattering and radiation forces, validated against analytical and semi-analytical solutions. The method accurately captures near- and far-field behavior across Rayleigh to geometrical optics regimes. Extensions to composite Janus particles reveal complex force, torque, and motion dynamics. Suitable for parallel simulations on complex geometries.

## 🧩 Related Publications

### [**Electromagnetic Scattering by Curved Surfaces and Calculation of Radiation Force: Lattice Boltzmann Simulations**](https://doi.org/10.1063/5.0234413)
*LBM applied to curved geometries for scattering and radiation-force computations.*

In this work, we apply the Lattice Boltzmann Method to electromagnetic scattering from curved surfaces and evaluate the resulting radiation forces. The method is systematically validated against analytical and semi-analytical solutions for conducting and dielectric cylinders, including corrugated and elliptical geometries. Excellent agreement is observed across different scattering regimes, establishing LBM as a reliable tool for accurate field and radiation-force computations on complex curved boundaries.

---

### [**Lattice Boltzmann Method for Electromagnetic Wave Scattering**](https://kwnsfk27.r.eu-west-1.awstrack.me/L0/https:%2F%2Fauthors.elsevier.com%2Fc%2F1mvnG_Wcs4XTA/1/0102019d76556fa0-4f3f5a5f-eacc-488a-b007-d61ea83cf099-000000/mCnP2aAsjGDLBiALlkiy8pPEhhQ=473)
*LBM as a numerical framework for broadband electromagnetic scattering.*

This paper proposes the Lattice Boltzmann Method as an alternative numerical framework for electromagnetic scattering computations. LBM is benchmarked across one-, two-, and three-dimensional configurations, covering plane-wave reflection and refraction, scattering from circular and hexagonal cylinders, and spherical geometries. The results are compared with analytical and Discretized Mie-Formalism solutions, demonstrating the accuracy, stability, and versatility of LBM for broadband electromagnetic simulations.

---

### [**Radiation Forces and Torques on Janus Cylinders**](https://arxiv.org/abs/2509.22308)
*LBM for composite materials: radiation force, torque, and motion of Janus cylinders.*

In this study, we extend the Lattice Boltzmann Method to model electromagnetic scattering and radiation-induced forces on composite (Janus) cylinders. LBM results for scattering and radiation forces/torques are benchmarked against semi-analytical solutions, showing excellent agreement. We further compute the coupled electromagnetic and hydrodynamic interactions to simulate the trajectories of Janus cylinders under combined radiation and viscous forces, highlighting their potential in optofluidic and active-matter systems.


