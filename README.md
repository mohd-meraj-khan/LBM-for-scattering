# LBM-for-scattering-of-EM-wave

## A 2D lattice Boltzmann solver for scattering and radiation force calculations

The solver utilizes the Lattice Boltzmann method (LBM) as described by [Hauser and Verhey](https://doi.org/10.1103/PhysRevE.96.063306). 
To optimize performance, the LBM code is implemented in C and interfaced using <b>ctypes</b> to execute within Python. Python was chosen 
due to its ease of coding and extensive library support.

Currently, the solver is designed to compute scattering width and radiation force for 2D geometries. It performs effectively across all three scattering regimes: Rayleigh, Mie, and geometrical optics.


### _Scattering width calculation_

Scattering width is a measure of how effectively an object scatters. For two dimensional scatterers it is defined as

$$\sigma = \lim_{r \to \infty} 2 \pi r \frac{|\pmb{ \mathcal{E}}^S|^2}{|\pmb{ \mathcal{E}}^I|^2}.$$

Here, $\sigma$ is the scattering width, $\pmb{ \mathcal{E}}^S$ and $\pmb{ \mathcal{E}}^I$ are the scattered and incident electric fields at a distance $r$ from the scatterer.



#### _Conducting cylinder_


#### _Dielectric cylinder_



### _Radiation force calculation_

Radiation force can be computed by integrating Maxwell's stress tensor over the surface of the scatterer.

$$\mathbb{T} = \varepsilon_0 \pmb{\mathcal{E}} \pmb{\mathcal{E}} + \mu_0 \pmb{\mathcal{H}} \pmb{\mathcal{H}} - \frac{1}{2} \left( \varepsilon_0 \pmb{\mathcal{E}} \cdot \pmb{\mathcal{E}} + \mu_0 \pmb{\mathcal{H}} \cdot \pmb{\mathcal{H}} \right) \mathbb{I}$$

Here, $\mathbb{T}$ is Maxwell's stress tensor, and $\mathbb{I}$ is the identity tensor.

The radiation force per unit length on the cylinder averaged over a time period can be computed by integrating average Maxwell's stress tensor over the perimeter of the cylinder
$$\langle {\bf F} \rangle = \oint \langle \mathbb{T} \rangle \cdot  \hat{{\bf n}} dl.$$

Here, $\hat{{\bf n}}$ is the unit normal vector at the perimeter of the cylinder.



#### _Conducting cylinder_

##### _Smooth circular_

##### _Corrugated elliptic_


#### _Dielectric cylinder_
