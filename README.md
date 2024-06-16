# LBM-for-scattering-of-EM-wave
A 2D lattice Boltzmann solver for scattering and radiation force calculations

<p>We usd the Lattice Boltzmann method given by <a href="https://doi.org/10.1103/PhysRevE.96.063306" target="_blank">Winners</a> </p>


<h2><i>Scattering width calculation</i></h2>

Scattering width is a measure of how effectively an object scatters. For two dimensional scatterers it is defined as

$$\sigma = \lim_{r \to \infty} 2 \pi r \frac{|\pmb{ \mathcal{E}}^S|^2}{|\pmb{ \mathcal{E}}^I|^2}.$$

Here, $\sigma$ is the scattering width, $\pmb{ \mathcal{E}}^S$ and $\pmb{ \mathcal{E}}^I$ are the scattered and incident electric fields at a distance $r$ from the scatterer.








<h2><i>Radiation force calculation</i></h2>



Radiation force can be computed by integrating Maxwell's stress tensor over the surface of the scatterer.

$$\mathbb{T} = \varepsilon_0 \pmb{\mathcal{E}} \pmb{\mathcal{E}} + \mu_0 \pmb{\mathcal{H}} \pmb{\mathcal{H}} - \frac{1}{2} \left( \varepsilon_0 \pmb{\mathcal{E}} \cdot \pmb{\mathcal{E}} + \mu_0 \pmb{\mathcal{H}} \cdot \pmb{\mathcal{H}} \right) \mathbb{I}$$

Here, $\mathbb{T}$ is Maxwell's stress tensor, and $\mathbb{I}$ is the identity tensor.



The radiation force per unit length on the cylinder averaged over a time period can be computed by integrating average Maxwell's stress tensor over the perimeter of the cylinder

$$\langle {\bf F} \rangle = \oint \langle \mathbb{T} \rangle \cdot  \hat{{\bf n}} dl.$$



Here, $\hat{{\bf n}}$ is the unit normal vector at the perimeter of the cylinder.
