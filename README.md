# LBM for scattering of EM wave

## A 2D lattice Boltzmann solver for scattering and radiation force calculations

The solver utilizes the Lattice Boltzmann method (LBM) as described by [Hauser and Verhey](https://doi.org/10.1103/PhysRevE.96.063306). 
To optimize performance, the LBM code is implemented in C and interfaced using <b>ctypes</b> to execute within Python. Python was chosen 
due to its ease of coding and extensive library support.

Currently, the solver is designed to compute scattering width and radiation force for 2D geometries. It performs effectively across all three scattering regimes: Rayleigh, Mie, and geometrical optics.


We compare LBM solutions with analytical solutions for smooth circular conducting and dielectric cylinders in scattering width and radiation force calculations. Additionally, we compare LBM solutions with [semi-analytical](https://doi.org/10.1364/OSAC.2.000298) solutions for corrugated elliptical conducting cylinders in radiation force calculations. In all cases, we find strong agreement between LBM solutions and both analytical and semi-analytical solutions across all three scattering regimes.


### _Scattering width calculation_

Scattering width is a measure of how effectively an object scatters. For two dimensional scatterers it is defined as

$$\sigma = \lim_{r \to \infty} 2 \pi r \frac{|\pmb{ \mathcal{E}}^S|^2}{|\pmb{ \mathcal{E}}^I|^2}.$$

Here, $\sigma$ is the scattering width, $\pmb{ \mathcal{E}}^S$ and $\pmb{ \mathcal{E}}^I$ are the scattered and incident electric fields at a distance $r$ from the scatterer.



#### _Conducting cylinder (smooth circular)_


The left side of Fig. 1 depicts the scattered electric field from the smooth circular conducting cylinder, while the right side shows bistatic (outer) and monostatic (inner) scattering widths for various $a / \lambda$ values across all three scattering regimes. In the right figure, solid lines represent analytical solutions, while dashed lines and markers denote LBM solutions.


<p align="center">
    <img src="https://github.com/mohd-meraj-khan/LBM-for-scattering/blob/main/Media/PEC_Es_2D_2_column.svg" style="display: inline-block;">
    <img src="https://github.com/mohd-meraj-khan/LBM-for-scattering/raw/main/Media/RCS_PEC.svg" style="display: inline-block;">
  <br>
  <em>Figure 1: Scattred electric field on left and bistatic (outset) and monostatic (inset) scattering width on right of the smooth circular conducting cylinder.
      Solid lines represents the analytical solutions whereas dashed lines and markers represents the LBM solutions.</em>
</p>







#### _Dielectric cylinder (smooth circular)_


We analyze smooth circular dielectric cylinders in two scenarios: (1) fixing the dielectric constant $\varepsilon_r$​ and varying the $a / \lambda$ values, and (2) fixing the $a / \lambda$ value and varying the dielectric constant $\varepsilon_r$​​ values.
 ### $\varepsilon_r = 2$


Here, we maintain a constant dielectric constant $\varepsilon_r$​ and vary the $a / \lambda$ values to encompass all three scattering regimes. The $a / \lambda$ values used are identical to those shown in Fig. 1


<p align="center">
    <img src="https://github.com/mohd-meraj-khan/LBM-for-scattering/blob/main/Media/RCS_Rayleigh_er_2.0.svg" style="display: inline-block;">
    <img src="https://github.com/mohd-meraj-khan/LBM-for-scattering/blob/main/Media/RCS_Mie_er_2.0.svg" style="display: inline-block;">
    <img src="https://github.com/mohd-meraj-khan/LBM-for-scattering/blob/main/Media/RCS_GO_er_2.0.svg" style="display: inline-block;">
  <br>
  <em>Figure 2: Bistatic (outset) and monostatic (inset) scattering width of the smooth circular dielectric cylinder for the Rayleugh, Mie and geometric optics regimes.
  Same legends as in Fig. 1.</em>
</p>

### $a/ \lambda = 0.5$


Here, we keep the $a / \lambda$ value fixed and vary the dielectric constant. In the case of a perfect conductor, the dielectric constant is infinite. As we increase the dielectric constant, the solution approaches that of a conducting cylinder.


<p align="center">
    <img src="https://github.com/mohd-meraj-khan/LBM-for-scattering/blob/main/Media/BRCS_ratio_0.5.svg" style="display: inline-block;">
    <img src="https://github.com/mohd-meraj-khan/LBM-for-scattering/blob/main/Media/MRCS_ratio_0.5.svg" style="display: inline-block;">
  <br>
  <em>Figure 3: Bistatic and monostatic scattering width of the smooth circular dielectric cylinder.</em>
</p>



### _Radiation force calculation_

Radiation force can be computed by integrating Maxwell's stress tensor over the surface of the scatterer.

$$\mathbb{T} = \varepsilon_0 \pmb{\mathcal{E}} \pmb{\mathcal{E}} + \mu_0 \pmb{\mathcal{H}} \pmb{\mathcal{H}} - \frac{1}{2} \left( \varepsilon_0 \pmb{\mathcal{E}} \cdot \pmb{\mathcal{E}} + \mu_0 \pmb{\mathcal{H}} \cdot \pmb{\mathcal{H}} \right) \mathbb{I}$$

Here, $\mathbb{T}$ is Maxwell's stress tensor, and $\mathbb{I}$ is the identity tensor.

The radiation force per unit length on the cylinder averaged over a time period can be computed by integrating average Maxwell's stress tensor over the perimeter of the cylinder
$$\langle {\bf F} \rangle = \oint \langle \mathbb{T} \rangle \cdot  \hat{{\bf n}} dl.$$

Here, $\hat{{\bf n}}$ is the unit normal vector at the perimeter of the cylinder.



#### _Conducting cylinder_



##### _Smooth circular_

Here, we present a comparison of the radiation force on smooth circular conducting cylinders as a function of $a / \lambda$. We consider $a / \lambda$ values ranging from 0.02 to 4, encompassing all three scattering regimes. The main figure illustrates the Mie and geometrical optics regimes, while the inset represents the Rayleigh regime.


<p align="center">
    <img src="https://github.com/mohd-meraj-khan/LBM-for-scattering/blob/main/Media/FxPEC.svg" style="display: inline-block;">
  <br>
  <em>Figure 4: Normalized radiation force on smooth circular conducting cylinder, Rayleigh regime (inset) and Mie and geometric optics regimes (outset).</em>
</p>


##### _Corrugated elliptic_


To assess the suitability of LBM for analyzing complex geometries, we refer to the work by [F. G. Mitri](https://doi.org/10.1364/OSAC.2.000298) , where the author employed the modal expansion method to compute the radiation force per unit length acting on corrugated and smooth elliptical conducting cylinders. Here, we have reproduce the radiation force calculations for corrugated cylinders using LBM.

Fig. 5 shows the total (incident + scattered) magnetic field for the corrugated elliptical conducting cylinders of different aspect ratios $'A'$.


<p align="center">
    <img src="https://github.com/mohd-meraj-khan/LBM-for-scattering/blob/main/Media/corrugatedSnap.svg" style="display: inline-block;">
  <br>
  <em>Figure 5: Total magnetic field due to corrugated conducting elliptical cylinders of different aspect ratios.</em>
</p>

Fig. 6 shows the comparison of the LBM solution with the semi-analytical solutions of the radiation force for the corrugated elliptical conducting cylinders shown in Fig. 5.

<p align="center">
    <img src="https://github.com/mohd-meraj-khan/LBM-for-scattering/blob/main/Media/corrugatedPEC.svg" style="display: inline-block;">
  <br>
  <em>Figure 6: Normalized radiation force for the cylinders shown in Fig. 5. Solid lines represent semi-analytical solution whereas the markers represents the LBM solutions</em>
</p>


#### _Dielectric cylinder (smooth circular)_

Here, we analyze the radiation force on smooth circular dielectric cylinders as a function of $a / \lambda$ for various dielectric constants, covering all three scattering regimes from $a / \lambda = 0.02$ to $a / \lambda = 2.6$. For clarity we have taken small window for each scattering regimes, the first three figures of Fig. 7. 



<p align="center">
    <img src="https://github.com/mohd-meraj-khan/LBM-for-scattering/blob/main/Media/FxRayleigh.svg" style="display: inline-block;">
    <img src="https://github.com/mohd-meraj-khan/LBM-for-scattering/blob/main/Media/FxMie.svg" style="display: inline-block;">
    <img src="https://github.com/mohd-meraj-khan/LBM-for-scattering/blob/main/Media/FxGO.svg" style="display: inline-block;">
    <img src="https://github.com/mohd-meraj-khan/LBM-for-scattering/blob/main/Media/FxExactLog.svg" style="display: inline-block;">
  <br>
  <em>Figure 7: Normalized radiation force on smooth circular dielectric cylinder.</em>
</p>


The sharp peaks in the above figure are called the resonance peaks where the amplitude of the fields are very high. Below we see three videos of total electric fields, where $a / \lambda = 0.95$ correspods to resonance peak. For all three videos, all the conditions are same except the value of $a / \lambda$. For all three cases, the amplitude of incident wave is 1. 

The sharp peaks observed in the above figure are resonance peaks, characterized by significantly high field amplitudes. Below, we present three videos displaying the total electric fields, with $a / \lambda = 0.95$ corresponding to a resonance peak. All conditions in three videos are identical, except for the value of $a / \lambda$, which varies across the three cases. In each video, the incident wave has an amplitude of 1.


### $$\varepsilon_r = 4, a/ \lambda = 0.91$$


https://github.com/mohd-meraj-khan/LBM-for-scattering/assets/153921085/97e00a97-a9b5-4946-8ff9-95362c6d3ef5



### $$\varepsilon_r = 4, a/ \lambda = 0.95$$


https://github.com/mohd-meraj-khan/LBM-for-scattering/assets/153921085/50ba4aed-5b85-46cf-b4e4-5827ca65818f



### $$\varepsilon_r = 4, a/ \lambda = 1.0$$


https://github.com/mohd-meraj-khan/LBM-for-scattering/assets/153921085/56dc300e-422f-4ac4-b635-2bb6ea705f3e







