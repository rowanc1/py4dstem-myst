---
title: Automated Crystal Orientation Mapping in py4DSTEM using Sparse Correlation Matching
short_title: ACOM Correlation
authors:
  - name: Colin Ophus
    corresponding: true
    email: cophus@gmail.com
    affiliations:
      - National Center for Electron Microscopy, Molecular Foundry, Lawrence Berkeley National Laboratory, 1 Cyclotron Road, Berkeley, CA, USA, 94720
  - name: Steven E Zeltmann
    affiliations:
      - Department of Materials Science and Engineering, University of California, Berkeley, CA,  94720
  - name: Alexandra Bruefach
    affiliations:
      - Department of Materials Science and Engineering, University of California, Berkeley, CA,  94720
  - name: Alexander Rakowski
    affiliations:
      - National Center for Electron Microscopy, Molecular Foundry, Lawrence Berkeley National Laboratory, 1 Cyclotron Road, Berkeley, CA, USA, 94720
  - name: Benjamin H Savitzky
    affiliations:
      - National Center for Electron Microscopy, Molecular Foundry, Lawrence Berkeley National Laboratory, 1 Cyclotron Road, Berkeley, CA, USA, 94720
  - name: Andrew M Minor
    affiliations:
      - National Center for Electron Microscopy, Molecular Foundry, Lawrence Berkeley National Laboratory, 1 Cyclotron Road, Berkeley, CA, USA, 94720
      - Department of Materials Science and Engineering, University of California, Berkeley, CA,  94720
  - name: MC Scott
    affiliations:
      - National Center for Electron Microscopy, Molecular Foundry, Lawrence Berkeley National Laboratory, 1 Cyclotron Road, Berkeley, CA, USA, 94720
      - Department of Materials Science and Engineering, University of California, Berkeley, CA,  94720
subject: Materials Science
doi: 10.48550/arXiv.2111.00171
arxiv: https://arxiv.org/abs/2111.00171
github: https://github.com/py4dstem/py4dstem_tutorials
license: CC-BY-4.0
math:
  '\ii': '{i\mkern1mu}'
  '\ft': '\mathcal{F}'
  '\ift': '\mathcal{F}^{-1}'
  '\tensor': '\overline{\overline{#1}}'
  '\pyFDSTEM': '\texttt{py4DSTEM}' # TODO: Link this https://github.com/py4dstem/py4DSTEM
  '\prismatic': '\texttt{Prismatic}' # TODO: Link this https://prism-em.com/
ext:
  si:
    Molar: M
---

+++ { "abstract": true }

## Abstract

Crystalline materials used in technological applications are often complex assemblies composed of multiple phases and differently oriented grains. Robust identification of the phases and orientation relationships from these samples is crucial, but the information extracted from the diffraction condition probed by an electron beam is often incomplete. We therefore have developed an automated crystal orientation mapping (ACOM) procedure which uses a converged electron probe to collect diffraction patterns from multiple locations across a complex sample. We provide an algorithm to determine the orientation of each diffraction pattern based on a fast sparse correlation method. We test the speed and accuracy of our method by indexing diffraction patterns generated using both kinematical and dynamical simulations. We have also measured orientation maps from an experimental dataset consisting of a complex polycrystalline twisted helical AuAgPd nanowire. From these maps we identify twin planes between adjacent grains, which may be responsible for the twisted helical structure. All of our methods are made freely available as open source code, including tutorials which can be adapted to perform ACOM measurements on diffraction pattern datasets.

+++

## Introduction

Polycrystalline materials are ubiquitous in technological applications. An ideal crystal structure can be fully defined with a small number of parameters: the 3 vectors defining its unit cell, and the position and species of each atom inside the unit cell {cite:p}`borchardt2011crystallography`. To fully describe crystalline materials in the real world however, we require a description of both the crystal lattice, and all defects present in a given material. These include point defects such as dopants, vacancies, or interstitials {cite:p}`dederichs1978lattice`, line defects such as dislocations {cite:p}`lesar2014simulations`, planar defects including internal boundaries and surfaces {cite:p}`tang2006diffuse`, and finally volume defects such as local strain fields {cite:p}`janssen2007stress`. One large subset of crystalline materials are polycrystalline phases, which consist of many small crystalline grains, arranged in either a random or organized fashion. Many material properties such as mechanical strength {cite:p}`thompson2000structure`, optical response {cite:p}`park2019efficient, londono2021intrinsic`, or thermal or electrical conductivity {cite:p}`castro2019role` are strongly modulated by the density and orientation of the boundaries between crystalline grains {cite:t}`thompson1995texture`. Thus characterizing the orientation of polycrystalline grains is essential to understanding these materials.

The two primary tools used to study the orientation of polycrystalline materials are electron backscatter diffraction (EBSD) in scanning electron microscopy (SEM), and transmission electron microscopy (TEM). EBSD can measure the orientation of crystalline grains with very high accuracy, but has limited resolution and is primarily sensitive to the surface of materials {cite:p}`humphreys2001review, wright2011review, wright2015introduction`. Alternatively, we can directly measure the atomic-scale structure and therefore the orientation of polycrystalline grains, either by using plane wave imaging in TEM {cite:p}`li2020constrained`, or by focusing the probe down to sub-atomic dimensions and scanning over the sample surface in scanning TEM (STEM) {cite:p}`peter2018segregation`. This is possible due to the widespread deployment of aberration correction for both TEM and STEM instruments {cite:p}`linck2016chromatic, ramasse2017twenty`. Atomic resolution imaging, however, strictly limits the achievable field-of-view, and requires relatively thin samples, and thus is primarily suited for measuring polycrystalline grain orientations of 2D materials {cite:p}`ophus2015large, qi2020near`.

Another approach to orientation mapping in TEM is to use diffraction space measurements. For crystalline materials, diffraction patterns will contain Bragg spots with spacing inversely proportional to the spacing of atomic planes which are approximately perpendicular to the beam direction (described by both the Laue condition and Bragg equations {cite:p}`fultz2012transmission`). To generate a spatially-resolved orientation map, we can focus a STEM probe down to dimensions of 0.5 to 50 nm, scan it over the sample surface, and record the diffraction pattern for each probe position, a technique referred to as nanobeam electron diffraction (NBED) {cite:p}`ozdol2015strain`, scanning electron nanobeam diffraction (SEND) {cite:p}`tao2009direct`, or four dimensional-scanning transmission electron microscopy (4D-STEM) (we choose this nomenclature for this text) due to the 4D shape of the collected data {cite:p}`bustillo20214d`. 4D-STEM experiments are increasingly enabled by fast direct electron detectors, as these cameras allow for much faster recording and much larger fields of view {cite:p}`ophus2019four, nord2020fast, paterson2020fast`.

By performing template matching of diffraction pattern libraries on 4D-STEM datasets, we can map the orientation of all crystalline grains with sufficient diffraction signal. This method is usually named automated crystal orientation mapping (ACOM), and has been used by many authors in materials science studies {cite:p}`zaefferer1994line, rauch2005rapid, kobler2013combination, maclaren2020comparison, londono20201d, jeong2021automated, zuo2021strategies`. ACOM experiments in 4D-STEM are highly flexible; two recent examples include {cite:t}`lang2021automated` implementing ACOM measurements in liquid cell experiments, and {cite:t}`wu2021seeing` adapting the ACOM method to a scanning confocal electron diffraction (SCED) experimental configuration. ACOM is also routinely combined with precession electron diffraction, where the STEM beam is continually rotated around a cone incident onto the sample, in order to excite more diffraction spots and thus produce more interpretable diffraction patterns {cite:p}`brunetti2011confirmation, moeck2011high, eggeman2015scanning`. Recently, {cite:t}`mehta2020unravelling` have combined simulations with machine learning segmentation to map orientations of 2D materials, and {cite:t}`yuan2021training` have used machine learning methods to improve the resolution and sensitivity of orientation maps by training on simulated data.

In this study, we introduce a new sparse correlation framework for fast calculation of orientation maps from 4D-STEM datasets. Our method is based on template matching of diffraction patterns along only the populated radial bands of a reference crystal's reciprocal lattice, and uses direct sampling of the first two Euler angles, and a fast Fourier transform correlation step to solve for the final Euler angle. We test our method on both kinematical calculations, and simulated diffraction experiments incorporating dynamical diffraction. Finally, we generate orientation maps of polycrystalline AuAgPd helically twisted nanowires, and use clustering to segment the polycrystalline structure, and map the shared (111) twin planes of adjacent grains.

## Methods

:::{figure} figure_ACOM_examples_v05.pdf
:name: Fig:ACOM_corr_match

**ACOM using correlation matching in $\pyFDSTEM$**. (a) Structure of fcc Au. (b) Atomic scattering factor of Au. (c) Structure factors for fcc Au. (d) Zone axes included in orientation plan. (e) Diffraction patterns for various orientations, and (f) corresponding orientation plan slices. (g) Correlogram maxima for each pattern in (e) as a function of zone axis, and (h) corresponding in-plane rotation correlation. Highest correlation scores are shown in (g) and (h) using red circles.
:::

### Structure Factor Calculations

The _structure factors_ of a given crystalline material are defined as the complex coefficients of the Fourier transform of an infinite crystal {cite:p}`spence1993accurate`. We require these coefficients in order to simulate kinematical diffraction patterns, and thus we briefly outline their calculation procedure here.

First, we define the reference crystal structure. This structure consists of two components, the first being its unit cell defined by its lattice vectors $\bm{a}$, $\bm{b}$, and $\bm{c}$ composed of positions in $\bm{r}=(x,y,z)$, the 3D real space coordinate system. The second component of a crystal structure is an array with dimensions $[N, 4]$ containing the fractional atomic positions $\bm{p}_n = (p_{\bm{a}}, p_{\bm{b}}, p_{\bm{c}})_n$ and atomic number $Z_n$, for the $n$th index of $N$ total atoms in the unit cell. Together these positions and atomic numbers are referred to as the atomic basis. Because $\bm{p}_n$ is given in terms of the lattice vectors, all fractional positions have values inside the range $[0, 1)$. The unit cell and real space Cartesian coordinates of the fcc Au structure are plotted in {numref}`Figure %s <Fig:ACOM_corr_match>`a.

All subsequent calculations are performed in reciprocal space (also known as Fourier space or diffraction space). Thus the next step is to compute the reciprocal lattice vectors, defined by {cite:p}`gibbs1884elements`

```{math}
\begin{eqnarray}
    \bm{a}^*
&=&
    \frac{\bm{b} \times \bm{c}}{
    \bm{a} \cdot [\bm{b} \times \bm{c}]}
    =
    \frac{\bm{b} \times \bm{c}}{V}
    \nonumber
    \\
    \bm{b}^*
&=&
    \frac{\bm{c} \times \bm{a}}{
    \bm{b} \cdot [\bm{c} \times \bm{a}]}
    =
    \frac{\bm{c} \times \bm{a}}{V}
    \nonumber
    \\
    \bm{c}^*
&=&
    \frac{\bm{a} \times \bm{b}}{
    \bm{c} \cdot [\bm{a} \times \bm{b}]}
    =
    \frac{\bm{a} \times \bm{b}}{V},
\end{eqnarray}
```

where $\times$ represents the vector cross product and $V$ is the cell volume in real space. Note that this definition does not include factors of $2\pi$, and therefore all reciprocal coordinates have spatial frequency units.

Next, we calculate the position of all reciprocal lattice points required for our kinematical diffraction calculation, given by

```{math}
\bm{g}_{hkl}
=
h \, \bm{a}^* +
k \, \bm{b}^* +
l \, \bm{c}^*,
```

where $h$, $k$, and $l$ are integers representing the reciprocal lattice index points corresponding to the Miller indices $(h,k,l)$. We include only points where $|\bm{q}_{hkl}| < q_{\rm{max}}$, where $\bm{q}=(q_x,q_y,q_z)$ are the 3D coordinates in reciprocal space, i.e.\ those which fall inside a sphere given by the maximum scattering vector $q_{\rm{max}}$. To find all reciprocal lattice coordinates, we first determine the shortest vector given by linear combinations of $(\bm{a}^*,\bm{b}^*,\bm{c}^*)$, and divide $q_{\rm{max}}$ by this vector length to give the range for $(h,k,l)$. We then tile $(h,k,l)$ in both the positive and negative directions up to this value, and then remove all points with vector lengths larger than $q_{\rm{max}}$.

The reciprocal lattice defined above represents all possible coordinates where the structure factor coefficients $V_g(\bm{q})$ could be non-zero. The structure factor coefficients depend only the atomic basis and are given by

```{math}
F_{hkl} =
\sum_{n=1}^N
f_n(|\bm{g}_{hkl}|)
\exp\left[
-2 \pi \ii (h,k,l) \cdot \bm{p}_n
\right],
```

where $f_n$ are the the single-atom scattering factors for the $n$th atom, which describe the scattering amplitude for a single atom isolated in space. There are multiple ways to parameterize $f_n$, but here we have chosen to use the factors defined by {cite:t}`lobato2014accurate` which are implemented in $\pyFDSTEM$. {numref}`Figure %s <Fig:ACOM_corr_match>`b shows the atomic scattering factor for an Au atom.

We have now defined all structure factor coefficients for a perfect infinite crystal as

```{math}
  V_g(\bm{q}) =
  \begin{cases}
    F_{hkl} & \text{if $\bm{q} = \bm{g}_{hkl}$}
    \\
    0 & \text{otherwise}.
  \end{cases}
```

{numref}`Figure %s <Fig:ACOM_corr_match>`c shows the structure factors of fcc Au, where the marker size denotes the intensity (magnitude squared) of the $F_{hkl}$ values.

### Calculation of Kinematical Diffraction Patterns

Here we briefly review the theory of kinematical diffraction of finite crystals, following {cite:t}`de2003introduction`. We can fully describe an electron plane wave by its wavevector $\bm{k}$, which points in the direction of the electron beam and has a length given by $|\bm{k}| = 1/\lambda$, where $\lambda$ is the (relativistically-corrected) electron wavelength. Bragg diffraction of the electron wave along a direction $\bm{k}'$ occurs when electrons scatter from equally spaced planes in the crystal, described in reciprocal space as

```{math}
:label: eq:diffraction_condition
\bm{k}' = \bm{k} + \bm{g}_{hkl}.
```

For elastic scattering, $\bm{k}'$ has the same length as $\bm{k}$, and so scattering can only occur along the spherical surface known as the _Ewald sphere construction_ {cite:p}`ewald1921berechunung`. This expression will almost never be satisfied by a perfect infinite crystal. However, real samples have finite dimensions, and thus the Fourier transform of their lattice will include a _shape factor_ $D(\bm{q})$ convolved with each reciprocal lattice point. Thus diffraction can still occur, as long as {eq}`eq:diffraction_condition` is approximately satisfied.

If the sample foil is tilted an angle $\alpha$ away from the beam direction, the vector between a reciprocal lattice point $\bm{g}$ and its closest point on the Ewald sphere has a length equal to

```{math}
:label: eq:excitation error
s_{\bm{g}} =
\frac{-\bm{g} \cdot(2 \bm{k} + \bm{g}) }
{2 | \bm{k} + \bm{g} | \cos(\alpha) }.
```

The $s_{\bm{g}}$ term is known as the _excitation error_ of a given reciprocal lattice point $\bm{g}$. When the excitation error $s_{\bm{g}}=0$, the Bragg condition is exactly satisfied. When the length of $s_{\bm{g}}$ is on the same scale as the extent of the shape factor, the Bragg condition is approximately satisfied.

A typical TEM sample can be approximately described as a slab or foil which is infinite in two dimensions, and with some thickness $t$ along the normal direction. The shape function of such a sample is equal to

```{math}
:label: Eq:shape_foil
D(q_z) = \frac{\sin(\pi q_z t)}{\pi q_z}.
```

Because this expression is convolved with each reciprocal lattice point, we can replace $q_z$ with the distance between the Ewald sphere and the reciprocal lattice point. For the orientation mapping application considered in this paper, we assume that $\alpha = 0$, and that the sample thickness $t$ is unknown. Instead, we replace {eq}`Eq:shape_foil` with the approximation

```{math}
:label: Eq:shape_gaussian
D(q_z) =
\exp\left(
-\frac{{q_z}^2}{2 \sigma^2}
\right),
```

where $\sigma$ represents the excitation error tolerance for a given diffraction spot to be included. We chose this expression for the shape function because it decreases monotonically with increasing distance between the diffraction spot and the Ewald sphere $q_z$, and produces smooth output correlograms.

To calculate a kinematic diffraction pattern for a given orientation $\bm{w}$, we loop through all reciprocal lattice points and use {eq}`eq:excitation error` to calculate the excitation errors. The intensity of each diffraction spot is given by the intensity of the structure factor $|F_{hkl}|^2$, reduced by a factor defined by either {eq}`Eq:shape_foil` or {eq}`Eq:shape_gaussian`. We define the position of the diffraction spots in the imaging plane by finding two vectors perpendicular to the beam direction, and projecting the diffraction vectors $q$ into this plane. The result is the intensity of each spot $I_m$, and its two spatial coordinates $(q_{m_x},q_{m_y})$, or alternatively their polar coordinates $q_m = \sqrt{{q_{m_x}}^2 + {q_{m_y}}^2}$ and $\gamma_m =$ arctan2$(q_{m_y}, q_{m_x})$. Note that the in-plane rotation angle is arbitrarily defined for kinematical calculations in the forward direction. The resulting diffraction patterns are defined by the list of $M$ Bragg peaks $(q_{m_x},q_{m_y},I_m)$ or $(q_m,\gamma_m ,I_m)$.

{numref}`Figure %s <Fig:ACOM_corr_match>`e shows diffraction patterns for fcc Au, along five different zone axes (orientation directions). Each pattern includes Bragg spots out to a maximum scattering angle of $q_{\rm{max}} = 1.5 \, \rm{\AA}^{-1}$, and each spot is labeled by the $(hkl)$ indices. The marker size shown for each spot scales with the amplitude of each spot's structure factor, decreased by {eq}`Eq:shape_gaussian` using $\sigma = 0.02 \, \rm{\AA}^{-1}$.

### Generation of an Orientation Plan

The problem we are solving is to identify the relative orientation between a given diffraction pattern measurement and a parent reference crystal. This orientation can be uniquely defined by a $[3 \times 3]$-size matrix $\tensor{\bm{m}}$, which rotates vectors $\bm{d}_0$ in the sample coordinate system to vectors $\bm{d}$ in the parent crystal coordinate system

```{math}
\begin{eqnarray}
    \begin{bmatrix}
       d_x \\
       d_y \\
       d_z
    \end{bmatrix}
    &=&
    \left[
      \begin{array}{ccc}
        u_x & v_x & w_x \\
        u_y & v_y & w_y \\
        u_z & v_z & w_z  \\
      \end{array}
    \right]
    \begin{bmatrix}
       d_{0_x} \\
       d_{0_y} \\
       d_{0_z}
    \end{bmatrix}
    \nonumber \\
    \bm{d}
    &=&
    \tensor{\bm{m}}
    \ \bm{d}_0,
\end{eqnarray}
```

where the first two columns of $\tensor{\bm{m}}$ given by $\bm{u}$ and $\bm{v}$ represent the orientation of the in-plane $x$ and $y$ axis directions of the parent crystal coordinate system, respectively, and the third column $\bm{w}$ defines the zone axis or out-plane-direction. The orientation matrix can be defined in many different ways, but we have chosen to use a $Z-X-Z$ Euler angle scheme {cite:p}`rowenhorst2015consistent`, defined as

```{math}
:label: Eq:rotation_matrices
\tensor{\bm{m}} =
\left[
  \begin{array}{ccc}
    C_1 & -S_1 & 0 \\
    S_1 & C_1 & 0 \\
    0 & 0 & 1  \\
  \end{array}
\right]
\left[
  \begin{array}{ccc}
    1 & 0 & 0  \\
    0 & C_2 & S_2 \\
    0 & -S_2 & C_2 \\
  \end{array}
\right]
\left[
  \begin{array}{ccc}
    C_3 & -S_3 & 0 \\
    S_3 & C_3 & 0 \\
    0 & 0 & 1  \\
  \end{array}
\right],
```

where $C_1 = \cos(\phi_1)$, $S_1 = \sin(\phi_1)$, $C_2 = \cos(\theta_2)$, $S_2 = \sin(\theta_2)$, $C_3 = \cos(\phi_3)$, and $S_3 = \sin(\phi_3)$. The Euler angles $(\phi_1, \theta_2, \phi_3)$ chosen are fairly arbitrarily, as are the signs of rotation matrices given above. This convention was chosen to preserve internal consistency and to produce orientation matrices with sorted directions.

In order to determine the orientation $\tensor{\bm{m}}$ of a given diffraction pattern, we use a two-step procedure. The first step is to calculate an _orientation plan_ $P((\phi_1, \theta_2), \phi_3, q_s)$ for a given reference crystal. The second step, which is defined in the following section, is to generate a _correlogram_ from each reference crystal, from which we directly determine the correct orientation.

The first two Euler angles $\phi_1$ and $\theta_2$ represent points on the unit sphere which will become the zone axis of a given orientation. The first step in generating an orientation plan is to select 3 vectors delimiting the extrema of the unique, symmetry-reduced zone axes possible for a given crystal. {numref}`Figure %s <Fig:ACOM_corr_match>`d shows these boundary vectors for fcc Au, which are given by the directions $[001]$, $[011]$, and $[111]$. We next choose a sampling rate or angular step size, and generate a grid of zone axes to test. We define a 2D grid of vectors on the unit sphere which span the boundary vectors by using spherical linear interpolation (SLERP) formula defined by {cite:t}`shoemake1985animating`. These points with a step size of $2^\circ$ are shown in {numref}`Figure %s <Fig:ACOM_corr_match>`d. The rotation matrices which transform the zone axis vector (along the z axis) are given by the matrix inverse of the first two terms in {eq}`Eq:rotation_matrices`.

We then examine the vector lengths of all non-zero reciprocal lattice points $\bm{g}_{hkl}$ and find all unique spherical shell radii $q_s$. These radii will become the first dimension of our orientation correlogram, where each radius is assigned one index $s$. We loop through all included zone axes, and calculate a polar coordinate representation of the kinematical diffraction patterns.

For each zone axis, the first step to compute the plan is to rotate all structure factor coordinates by the matrix inverse of the first two terms in {eq}`Eq:rotation_matrices`. Next, we compute the excitation errors $s_g$ for all peaks assuming a $[0,0,1]$ projection direction, and the in-plane rotation angle of all peaks $\gamma_{\bm{q}}$. The intensity values of the orientation plan for all $q_s$ shells and in-plane rotation values $\phi_3$ are defined using the expression

```{math}
\begin{eqnarray}
    && P_0((\phi_1, \theta_2), \phi_3, q_s) =
    \sum_{ \left\{ \bm{g} \,:\, |\bm{g}| = q_s  \right\} }
    {q_s}^\gamma |V_{\bm{g}}|^\omega \times
    \nonumber \\
    &&
    \rm{max} \left\{
        1 - \frac{1}{\delta} \sqrt{
        s_{\bm{g}}^2 +
        \left[
        \rm{mod}(
        \phi_3 - \gamma_{\bm{g}} + \pi, 2 \pi) - \pi
        \right]^2 {q_s}^2
        }, 0
    \right\},
    \nonumber
\end{eqnarray}
```

where $\delta$ is the correlation kernel size, $\gamma$ and $\omega$ represent the power law scaling for the radial and peak amplitude terms respectively, $\rm{max}(...)$ is the maximum function, which returns the maximum of its two arguments, $\rm{mod}(...)$ is the modulo operator, and the summation includes only those peaks $\bm{g}$ which belong to a given radial value $q_s$. We have used the combined indexing notation for $(\phi_1, \theta_2)$ to indicate that in practice, this dimension of the correlation plan contains all zone axes, and thus the entire array has only 3 dimensions. The correlation kernel size $\delta$ defines the azimuthal extent of the correlation signal for each reciprocal lattice point. Note that {eq}`Eq:shape_gaussian` and {eq}`Eq:shape_foil` are not used for the calculation of orientation plans.

We normalize each zone axis projection using the function

```{math}
\nonumber
A(\phi_1, \theta_2)
= \frac{1}{\sqrt{
\sum_{\phi_3}
\sum_{q_s}
{P_0((\phi_1, \theta_2), \phi_3, q_s)}^2
}},
```

yielding the final normalized orientation plan

```{math}
P((\phi_1, \theta_2), \phi_3, q_s) = A(\phi_1, \theta_2) P_0((\phi_1, \theta_2), \phi_3, q_s)
```

By default, we have weighted each term in the orientation plan with the prefactor $q_s |V_g|$, i.e.\ setting $\gamma=\omega=1$. The $q_s$ term gives slightly more weight to higher scattering angles, while the $|V_g|$ term is used to weight the correlation in favour of peaks with higher structure factor amplitudes, which was found to be more reliable than weighting the orientation plan by $|V_g|^2$, which weights each peak by its structure factor intensity.

{numref}`Figure %s <Fig:ACOM_corr_match>`f shows 2D slices of the 3D orientation plan, for the 5 diffraction patterns shown in {numref}`Figure %s <Fig:ACOM_corr_match>`e. The in-plane rotational symmetry of each radial band is obvious for the low index zone axes, e.g. for the $[001]$ orientated crystal, the first row of the corresponding orientation plan consists of four spots which maintains the 4-fold symmetry of the diffraction pattern and can be indexed as $[020]$, $[200]$, $[0\overline{2}0]$ and $[\overline{2}00]$. The final step is to take the 1D Fourier transform along the $\phi_3$ axis in preparation for the Fourier correlation step defined in the next section.

### Correlation Pattern Matching

For each diffraction pattern measurement, we first measure the location and intensity of each Bragg disk by using the template matching procedure outlined by {cite:t}`savitzky2021py4dstem`. The result is a set of $M$ experimental diffraction peaks defined by $(q_m,\gamma_m ,I_m)$. From these peaks, we calculate the sparse polar diffraction image $X(\phi_3, q_s)$ using the expression

```{math}
:label: eq:image_polar
\begin{eqnarray}
    && X(\phi_3, q_s) =
    \sum_{\{ q_m \,: \, |q_m - q_s| < \delta \}}
    {q_m}^\gamma {I_m}^{\omega/2} \times
    \rm{max} \left\{
        1 - \right.
    \\ % \nonumber
    &&
    \left. \frac{1}{\delta} \sqrt{
    (q_m - q_s)^2 +
    \left[
    \rm{mod}(
    \phi_3 - \gamma_m + \pi, 2 \pi) - \pi
    \right]^2 {q_s}^2
    }, 0
    \right\}.
    \nonumber
\end{eqnarray}
```

By default, we again use prefactors weighted by the peak radius and estimated peak amplitude given by the square root of the measured disk intensities. However, if the dataset being analyzed contains a large number of different sample thicknesses, multiple scattering can cause strong oscillations in the peak amplitude values. As we will see in the simulations below, in these situations the best results may be achieved by setting $\omega=0$, i.e. ignoring peak intensity and weighting only by the peak radii. Note that in the diffraction image, the correlation kernel size $\delta$ again gives the azimuthal extent of the correlation signal. However, in {eq}`eq:image_polar` it also sets the range over which peaks are included in a given radial bin, and the fraction of the intensity assigned to each radial bin.

Next, we calculate the correlation $C((\phi_1, \theta_2), \phi_3)$ of this image with the orientation plan using the expression

```{math}
\begin{eqnarray}
    && C((\phi_1, \theta_2), \phi_3)
    =
    \sum_{q_s}
    \ift \left\{  \right.
    \nonumber \\
    &&
    \left.
        \ft\left\{
        {P((\phi_1, \theta_2), \phi_3, q_s)}
        \right\}^*
        \ft\left\{
        X(\phi_3, q_s)
        \right\}
    \right\},
    \nonumber
\end{eqnarray}
```

where $\ft$ and $\ift$ are 1D forward and inverse fast Fourier transforms (FFTs) respectively along the $\phi_3$ direction, and the ${}^*$ operator represents taking the complex conjugate. We use this correlation over $\phi_3$ to efficiently calculate the in-plane rotation of the diffraction patterns. The maximum value in the correlogram will ideally correspond to the most probable orientation of the crystal. In order to account for mirror symmetry of the 2D diffraction patterns, we can also compute the correlation

```{math}
\begin{eqnarray}
    && C_{\rm{mirror}}((\phi_1, \theta_2), \phi_3)
    =
    \sum_{q_s}
    \ift \left\{  \right.
    \nonumber \\
    &&
    \left.
        \ft\left\{
        {P((\phi_1, \theta_2), \phi_3, q_s)}
        \right\}^*
        \ft\left\{
        X(\phi_3, q_s)
        \right\}^*
    \right\},
    \nonumber
\end{eqnarray}
```

where the mirror operation is accomplished by taking the complex conjugate of $\ft\left\{X(\phi_3, q_s)\right\}$. For each zone axis $(\phi_1, \theta_2)$, we take the maximum value of $C$ and $C_{\rm{mirror}}$ in order to account for this symmetry. {numref}`Figure %sg and %sh <Fig:ACOM_corr_match>` show 5 output correlograms, for the 5 diffraction patterns shown in {numref}`Figure %s <Fig:ACOM_corr_match>`e. For each zone axis $(\phi_1, \theta_2)$, we have computed the maximum correlation value, which are plotted as a 2D array in {numref}`Figure %s <Fig:ACOM_corr_match>`g. In each case, the highest value corresponds to the correct orientation.

{numref}`Figure %s <Fig:ACOM_corr_match>`h shows the correlation values along the $\phi_3$ axis, for the $(\phi_1, \theta_2)$ bins with the highest correlation value in {numref}`Figure %s <Fig:ACOM_corr_match>`g. The symmetry of the correlation values in {numref}`Figure %s <Fig:ACOM_corr_match>`h reflect the symmetry of the underlying patterns. For the $[0,0,1]$, $[0,1,1]$, and $[1,1,1]$, diffraction patterns, the in-plane angle $\phi_3$ correlation signals have 4-fold, 2-fold, and 6-fold rotational symmetry respectively. By contrast, the asymmetric diffraction patterns with zone axes $[1,1,3]$, and $[1,3,5]$ have only a single best in-plane orientation match.

### Matching of Overlapping Diffraction Patterns

In order to match multiple overlapping crystal signals, we have implemented an iterative detection process. First, we use the above algorithm to determine the best fit orientation for a given pattern. Next, the forward diffraction pattern is calculated for this orientation. We then loop through all experimental peaks, and any within a user-specified deletion radius are removed from the pattern. Peaks which are outside of this radius, but within the correlation kernel size, have their intensities reduced by a factor defined by the linear distance between the experimental and simulated peaks divided by the distance between the correlation kernel size and the deletion radius. Then, the ACOM correlation matching procedure is repeated until the desired number of matches have been found, or no further orientations are found. Note that while we could update the correlation score after peak deletion, we output the original magnitude of the full pattern correlogram in order to accurately calculate the probability of multiple matches.

:::{figure} figure_orientation_plans_v04.pdf
:name: Fig:ACOM_plans

Examples of alternative orientation plan types in $\pyFDSTEM$. Fiber texture examples where (a) orientations fully orbit around a single zone axis (the fiber axis), or (b) contain only a symmetry-reduced wedge of zone axes which orbit around a the fiber axis. (c) Examples of orientation plans generated directly from Materials Project entries {cite:p}`jain2013commentary`, using `pymatgen` symmetries {cite:p}`ong2013python`.
:::

### ACOM Integration into $\pyFDSTEM$

The ACOM pattern matching described has been implemented into the $\pyFDSTEM$ python toolkit written by {cite:t}`savitzky2021py4dstem`. A typical ACOM workflow starts with using $\pyFDSTEM$ to import the 4D dataset and one or more images of the vacuum probe. We then use a correlation template matching procedure to find the positions of all diffracted disks at each probe position {cite:p}`pekin2017optimizing`. We use the correlation intensity of each detected peak as an estimate of the peak's intensity. The resulting set of $M$ peaks defined by the values $(q_m, \gamma_m, I_m)$ are stored as a _PointList_ object in $\pyFDSTEM$. Because the number of peaks detected at each probe position can vary, we store the full set of all detected peaks in a _PointListArray_ object in $\pyFDSTEM$, which provides an interface to the ragged structured numpy data.

Most experimental datasets contain some degree of ellipticity, and the absolute pixel size must be calibrated. We perform these corrections on the set of measured diffraction disks using the $\pyFDSTEM$ calibration routines defined by {cite:t}`savitzky2021py4dstem`. We know that the correlation approach is relatively robust against both ellipticity and small errors in the reciprocal space pixel size. However, precise phase mapping may require us to distinguish between crystals with similar lattice parameters; these experiments will require accurate calibration.

We perform ACOM in $\pyFDSTEM$ by first creating a _Crystal_ object, either by specifying the atomic basis directly, or by using the `pymatgen` package {cite:p}`ong2013python` to import structural data from crystallographic information files (CIF), or the Materials Project database {cite:p}`jain2013commentary`. The _Crystal_ object is used to calculate the structure factors, and generate an orientation plan. The final step is to use the orientation plan to determine the best match (or matches) for each probe position, from the list of calibrated diffraction peaks. If the sample contains multiple phases, we perform the orientation plan calculation and correlation matching for each unique crystal structure.

In addition to specifying the orientation plan spanning 3 vectors as in {numref}`Figure %s <Fig:ACOM_corr_match>`, we define additional methods to describe the space of possible orientations. One such example is _fiber texture_, where we assume the crystals are all orientated near a single zone axis known as the fiber axis, shown in {numref}`Figures a and b <Fig:ACOM_plans>`. We can vary the angular range of zone axes included away from the fiber axis as in {numref}`Figure %s <Fig:ACOM_plans>`a, as well as choose the azimuthal range around this axis as in {numref}`Figure %s <Fig:ACOM_plans>`b to account for symmetry around the fiber axis. Alternatively, an "automatic" option is provided, which uses `pymatgen` to determine the symmetry of the structure and automatically choose the span of symmetrically unique zone axes which should be included in the orientation plan, based on the point group symmetry {cite:t}`de2003introduction`. This is shown for a selection of different Materials Project database entries in {numref}`Figure %s <Fig:ACOM_plans>`c.

### Simulations of Diffraction Patterns from Thick Samples

One important metric for the performance of an orientation mapping algorithm is how well it performs when the diffraction patterns contain significant amounts of multiple scattering. We have therefore used our ACOM algorithm to measure the orientation of simulated diffraction patterns from samples tilted along many directions, over a wide range of thicknesses. We performed these simulations using the multislice algorithm {cite:p}`cowley1957scattering`, and methods defined by {cite:t}`kirkland2020advanced` and {cite:t}`ophus2017fast`. These methods are implemented in the `prismatic` simulation code by {cite:t}`dacosta2021prismatic`. The diffraction patterns were generated using a acceleration potential of 300 keV, a 0.5 mrad convergence angle, with real space and reciprocal pixel sizes of 0.05 $\rm{\AA}$ and 0.01 $\rm{\AA}^{-1}$ respectively, with 4 frozen phonons. In total we have simulated 3750 diffraction patterns from Cu, Ag, and Au fcc crystals, over 25 zone axes ($[0,0,1]$ to $[3,4,4]$ excluding symmetrically redundant reflections) and thicknesses up to 100 nm with a 2 nm step size. These diffraction patterns were generated using the simulation pipeline and database defined by {cite:t}`rakoski2021database`.

### Chemical Synthesis of Twisted AuAgPd Nanowires

The performed synthesis was modified from a known method given by {cite:t}`wang2011`. We prepared the following solutions: {si}`500 <\milli\Molar>` PVP (MW 40,000) in DMF, {si}`50 <\milli\Molar>` {chem}`HAuCl4` in DMF, {si}`50 <\milli\Molar>` {chem}`AgNO3` in MilliQ water, and {si}`400 <\milli\Molar>` L-ascorbic acid in MilliQ water. We created the reaction solution in a 4 mL vial (washed 3x with MilliQ water and acetone) by mixing {si}`800 <\micro\liter>` DMF, {si}`100 <\micro\liter>` PVP, {si}`20 <\micro\liter>` {chem}`HAuCl4`, and {si}`20 <\micro\liter>` {chem}`AgNO3`. We vortexed the solution, then added {si}`100 <\micro\liter>` of L-ascorbic acid solution drop-wise to the mixture while gently swirling. At this point, the color changed from pale yellow to clear. We left the solution at room temperature for 7 days, at which point the solution was light brown/purple. The primary product of this reaction was straight, ultrathin Au-Ag nanowires ({si}`2 <\nano\meter>` in diameter).

To twist the underlying ultrathin Au-Ag nanowires, we prepared solutions of {si}`1.875 <\milli\Molar>` L-ascorbic acid and {si}`2 <\milli\Molar>` {chem}`H2PdCl2` in MilliQ water. In a 4 mL vial (3x washed with MilliQ water/acetone), we added {si}`50 <\micro\liter>` of the Au-Ag reacted solution to {si}`640 <\micro\liter>` of the L-ascorbic acid solution. Finally, we added {si}`60 <\micro\liter>` of the {chem}`H2PdCl4` solution and allowed the sample to incubate for at least 30 minutes. We purified the reaction solution by centrifuging the product down at 7500 rpm for 4 minutes. We decanted the supernatant, and then rinsed the reaction with MilliQ water 3 times and re-dispersed in MilliQ water. We prepared TEM samples of this material by depositing {si}`10 <\micro\liter>` of purified nanowire solution onto 400 mesh formvar/ultrathin carbon grids.

### 4D-STEM Experiments with Patterned Apertures

We collected the experimental data using a double aberration-corrected modified FEI Titan 80-300 microscope (the TEAM I instrument at the National Center for Electron Microscopy within Lawrence Berkeley National Laboratory). This microscope is equipped with a Gatan K3 detector and Continuum spectrometer, and was set to collect diffraction patterns integrated over 0.05 seconds, with 4x binning. We used an accelerating voltage of 300 keV, an energy slit of 20 eV, and a spot size of 6. We used a {si}`10 <\micro\meter>` bullseye aperture to form the STEM probe in order to improve detection precision of the Bragg disks {cite:p}`zeltmann2020`. We used a convergence semiangle of {si}`2 <\milli\radian>`, with a camera length of 1.05 m. We recorded the experimental dataset using a step size of {si}`5 <\angstrom>`, with a total of 286 and 124 steps in the x and y directions.

:::{figure} figure_ACOM_Au_error_v03.pdf
:name: Fig:ACOM_kinematical
**Zone axis misorientation as a function of sampling and maximum scattering angle for kinematical simulations.** The mean tilt error and number of patterns matched per second are shown inset for each panel.
:::

:::{figure} figure_sim_03.pdf
:name: Fig:dynamical_diffractions

**Dynamical simulated diffraction patterns.** (a) Example diffraction patterns for Au oriented to the [011] zone axis for 10-80 nm thick slices. (b) Plots showing the mean zone axis misorientation in degrees as a function of thickness for Cu, Ag, and Au. Each plot shows the errors for correlation prefactors of $q_{\rm{s}} |V_{\rm{g}}|$ (red) and $q_{\rm{s}}$ (blue).
:::

## Results and Discussion

### ACOM of Kinematical Calculated Diffraction Patterns

For the first test of our correlation method, we applied it to the same patterns calculated to generate an orientation plan for fcc Au. Next, we measured the calculation time and angular error between the measured and ground truth zone axes for each pattern. These results are plotted in {numref}`Figure %s <Fig:ACOM_kinematical>` for 3 different maximum scattering angles $k_{\rm{max}}$, and angular sampling of $1^\circ$ and $2^\circ$.

The results in {numref}`Figure %s <Fig:ACOM_kinematical>` show that the angular error in zone axis orientation is relatively insensitive to the angular sampling. However, the angular error drops by a factor of 10 from approximately $\approx 3^\circ$ to $\approx 0.3^\circ$ when increasing the maximum scattering angle included from $k_{\rm{max}} = 1 \, {\rm{\AA}}^{-1}$ to $1.5 \, {\rm{\AA}}^{-1}$, and by another factor of 2-3 when increasing $k_{\rm{max}}$ to $2 \, {\rm{\AA}}^{-1}$. This is unsurprising, as examining {numref}`Figure %s <Fig:ACOM_corr_match>`e shows that there is a large increase in the number of visible Bragg spots outside of $k_{\rm{max}} = 1 \, {\rm{\AA}}^{-1}$, and because Bragg disks at higher scattering angles provide better angular precision relative to low $k$ disks. This result emphasizes the importance of collecting as wide of angular range as possible when performing orientation matching of 4D-STEM data.

The inset calculation times reported are for the single-threaded ACOM implementation in $\pyFDSTEM$, running in Anaconda {cite:p}`anaconda` on a laptop with an Intel Core i7-10875H processor, running at 2.30 GHz. The calculation times can be increased by an order of magnitude or more when running in parallel, or by using a GPU to perform the matrix multiplication and Fourier transform steps.

### ACOM of Dynamical Simulated Diffraction Patterns

In diffraction experiments using thick specimens, the electron beam can scatter multiple times, a phenomenon known as dynamical diffraction. This effect is especially pronounced in diffraction experiments along low index zone axes, where the diffracted peak intensities oscillate a function of thickness. In order to test the effect of oscillating peak intensities on our ACOM method, we have simulated diffraction patterns for Cu, Ag, and Au fcc crystals, along multiple zone axes. Some example diffraction patterns for the [011] zone axis of Au are plotted in {numref}`Figure %s <Fig:dynamical_diffractions>`a. We see that all diffraction spots have intensities which oscillate multiple times as a function of thickness.

We performed ACOM by generating orientation plans with an angular sampling of $2^\circ$, a correlation kernel size of 0.08 $\rm{\AA}^{-1}$, and maximum scattering angles of $k_{\rm{max}}$ = 1.0, 1.5, and 2.0 $\rm{\AA}^{-1}$. We kept the radial prefactor of weighting set to $\gamma $ = 1, and tested peak amplitude prefactors of $\omega$ = 1.0, 0.5 and 0.0. The average zone axis angular misorientation as a function of thickness is plotted in {numref}`Figure %s <Fig:dynamical_diffractions>`b. In total, we performed orientation matching on 3750 diffraction patterns, and a total of 33750 correlation matches on a workstation with an AMD Ryzen Threadripper 3960X CPU (2.2 GHz, baseclock). The typical number of patterns matched per second were of ~80-90, 45-55 and 25-30 patterns/s for $k_{\rm{max}}$ values of 1.0, 1.5, and 2.0 $\rm{\AA}^{-1}$ respectively.

As expected, the errors are higher than those achieved under kinematic conditions, and the trend for smaller errors with larger $k_{\rm{max}}$ is also preserved (mean errors of 7.25$^\circ$, 3.09$^\circ$ and 1.39$^\circ$ for $k_{\rm{max}}$ values 1.0, 1.5 and 2.0 $\rm{\AA}^{-1}$ respectively, $\gamma $ = 1, $\omega$ = 0.25). We did not observe any dependence of the orientation accuracy on the simulation thickness. Despite the correlation prefactor $|Vg|$ performing well for the examples shown in {numref}`Figure %s <Fig:ACOM_corr_match>`, for the dynamical diffraction simulations along zone axes it was out-performed by prefactors of both $\sqrt{|V_g|}$ ($\omega$ = 0.5) and omitting the peak amplitude prefactor altogether ($\omega$ = 0). We therefore suggest that when mapping samples with a large range of thicknesses, or many crystals aligned to low index zone axes, the position of the diffracted peaks is significantly more important than their amplitudes or intensities. One possible method to increase the accuracy while using higher amplitude prefactors is to perform an experiment which recovers more kinematical values for the diffracted peak intensities, for example by precessing the electron beam when recording diffraction patterns {cite:p}`midgley2015precession, jeong2021automated`. We note that there is likely no global optimal choice of orientation mapping hyperparameters for all materials and thicknesses, and this may be a worthwhile topic for future investigations.

### 4D-STEM ACOM of Twisted AuAgPd Nanowires

We have tested our ACOM algorithm with a 4D-STEM dataset collected for an AuAgPd nanowires. An image of the vacuum bullseye STEM probe is shown in {numref}`Figure %s <Fig:exp_structure>`a. For each detector pixel, we have calculated the maximum value across all STEM probe positions to generate a _maximum diffraction pattern_, shown in {numref}`Figure %s <Fig:exp_structure>`b. The beamstop used to block the center beam is visible, as well as various crystalline diffraction rings out to approximately 1.4 ${\rm{\AA}}^{-1}$.

:::{figure} figure_AuAg_exp_v01.pdf
:name: Fig:exp_structure
**4D-STEM scan of twisted polycrystalline AuAgPd nanowires.** (a) Diffraction image of probe over vacuum, showing bullseye pattern. (b) Maximum of each pixel in diffraction space over all probe positions. (c) Histogram of all peak locations detected by correlation in $\pyFDSTEM$ of (a) with each pattern included in (b). (d) HAADF-STEM image of the sample. (e) 1D histogram of scattering vectors, with fcc AuAg inverse plane spacings overlaid.
:::

After performing the correlation peak finding algorithm in $\pyFDSTEM$, we have an estimated position and intensity of all detected Bragg peaks. A 2D histogram of these peaks, known as a Bragg vector map, is plotted in {numref}`Figure %s <Fig:exp_structure>`c. Sharp polycrystalline diffraction rings are clearly visible, as well as false positives generated by the beamstop edge. These false positives were manually removed by using a mask generated from an image of the beamstop. A high angle annular dark field (HAADF) image was simultaneously recorded during the 4D-STEM data collection, which is shown in {numref}`Figure %s <Fig:exp_structure>`d.

The final experimental pre-processing steps are to calibrate the diffraction pattern center, the elliptical distortions, and the absolute pixel size. We performed these steps by fitting an ellipse to the $(022)$ diffraction ring, and by assuming a lattice constant of $4.08 \, \rm{\AA}$, corresponding to the fcc Au structure {cite:p}`maeland1964lattice`. This process is explained in more detail by {cite:t}`savitzky2021py4dstem`. We assumed that the Ag lattice constant is similar to that of Au. Despite the presence of Pd in the nanowires, there was no significant presence of secondary grains corresponding to the smaller lattice of fcc Pd grains. An intensity histogram of the corrected Bragg peak scattering angles are shown in {numref}`Figure %s <Fig:exp_structure>`e. We have overlaid the 5 smallest scattering angles of Au on {numref}`Figure %s <Fig:exp_structure>`e to show the accuracy of the correction.

:::{figure} figure_AuAg_orientation_v03.pdf
:name: Fig:ACOM_AuAgPd
**Orientation mapping of polycrystalline AuAgPd nanowires.** (a) Total of measured correlation signal for each probe position. (b) Estimated number of patterns indexed for each probe position. (c) Example of 2 orientations indexed from a single diffraction pattern, collected at the position indicated by the arrow shown in (b). (d) Orientation maps of the 3 highest correlation signals for each probe position. Legend shown above.
:::

We have performed ACOM on the AuAgPd nanowire sample, with the results shown in {numref}`Figure %s <Fig:ACOM_AuAgPd>` shown for up to 3 matches for each diffraction pattern. For each probe position, the sum of the maximum detected correlation signals for up to three matches are shown in {numref}`Figure %s <Fig:ACOM_AuAgPd>`a. The structure is in good agreement with {numref}`Figure %s <Fig:exp_structure>`, though with additional modulations due to some grains generating more diffraction signal than others. Using a correlation intensity threshold of 0.5, we have plotted the number of matching patterns in {numref}`Figure %s <Fig:ACOM_AuAgPd>`b. The threshold of 0.5 was arbitrary chosen as a lower bound for a potential match. Examples of 2 matches to a single diffraction pattern are plotted in {numref}`Figure %s <Fig:ACOM_AuAgPd>`c. In this figure the correlation score for the first matched pattern was higher than the second. The second match found shows some deformation between the measured and simulated Bragg peak positions, and matches less peaks. It therefore produces a lower correlation score, which can be used to threshold the results as in {numref}`Figure %s <Fig:ACOM_AuAgPd>`d.

{numref}`Figure %s <Fig:ACOM_AuAgPd>`d shows the 3D orientations for all probe positions, with the 3 best matches shown. Each image is masked by the total correlation signal, so that low correlation values are colored black. Almost every diffraction pattern with Bragg disks detected was indexed for at least one orientation with high confidence. Additionally, the patterns are very consistent, with a large number of adjacent probe positions recording the same orientation. Some secondary grains are also clearly visible in the second-best match, while very few patterns have been assigned a third match with high confidence.

:::{figure} figure_AuAg_twins_v08.pdf
:name: Fig:ACOM_111_twins
**Orientation analysis of grains in AuAgPd nanowires.** (a) Crystal grains, with in-plane (111) planes colored by orientation. (b) (111) planes shared by two overlapping grains.
:::

In order to investigate the grain organization of the AuAgPd nanowires, we have performed clustering analysis on the orientation maps. Grains with similar orientations have been clustered together by looping through each probe position and comparing its orientation to its neighbors. Grains with at least 10 contiguous probe positions are shown in {numref}`Figure %s <Fig:ACOM_111_twins>`a. $(111)$ planes which lie in the image plane are overlaid onto the grain strucure, colored by their orientation. Confirming our observations in {numref}`Figure %s <Fig:ACOM_AuAgPd>`d, only a few grains with substantial overlap were reliably identified. This might be due to the low thickness of the sample (only a single grain along the beam direction), some grains not being oriented close enough to a zone axis to be detected, or multiple scattering deviations in the diffracted signal. There is a noticable bias in the orientation of the $(111)$ planes, which tend to be oriented horizontally near the growth direction of the nanowires.

%There does not appear to be any noticeable pattern in the overall orientations of the grains.

One hypothesis for the growth mode of these twisted nanowires is that adjacent grains are connected by $(111)$ twin planes, forming local helical structures to give the observed twisted structures. To test this hypothesis, we determined the position of $(111)$ planes from {numref}`Figure %s <Fig:ACOM_111_twins>`a which are shared by two overlapping grains. {numref}`Figure %s <Fig:ACOM_111_twins>`b shows the location of these shared $(111)$ planes (with plane normal differences below $8^\circ$), colored by the normal vector of the plane. Many shared $(111)$ planes were detected, most with normal vectors aligned to the wire growth direction. These observations support the hypothesis that these nanowires are composed of grains connected by $(111)$ twin planes.

%This is reasonable since in order to detect these shared planes, we require the grains on either side to be oriented along a low enough index zone axis. If both grains have a $(111)$ plane normal in the sample plane, the odds of both being in a good diffraction condition are increased.

These experimental observations demonstrate the efficacy of our ACOM method. In order to improve these results, we will need to collect diffraction data with a wider angular range. This can be achieved by using precession electron diffraction {cite:p}`rouviere2013improved`, multibeam electron diffraction {cite:t}`hong2021multibeam`, or by tilting the sample or beam and recording multiple 4D-STEM datasets {cite:p}`meng2016three`.

## Conclusion

We have introduced an efficient and accurate method to perform automated crystal orientation mapping, using a sparse correlation matching procedure. We have implemented our methods into the open source $\pyFDSTEM$ toolkit, and demonstrated the accuracy of our method using simulated diffraction patterns. We also applied ACOM to an experimental scan of a complex helical polycrystalline nanowire, where we were able to identify shared twin planes between adjacent grains which may be responsible for the twisted helical geometry. All of our methods have been made freely available to the microscopy community as open source codes. We believe that our implementation of ACOM is efficient and accurate enough to be Incorporated into automated online TEM software {cite:p}`spurgeon2021towards`. In the future, we will improve our ACOM method using machine learning methods {cite:t}`munshi2021ml`, and we will extend our ACOM methods to include multibeam electron diffraction experiments {cite:p}`hong2021multibeam`.

## Source Code and Data Availability

All code used in this manuscript is available on the [py4DSTEM GitHub repository](https://github.com/py4dstem/py4DSTEM/tree/acom), and the tutorial notebooks are available on the [py4DSTEM tutorial repository](https://github.com/py4dstem/py4DSTEM_tutorials/tree/main/notebooks/acom). All simulated and experimental 4D-STEM datasets are available at [links will be added after publication].

## Acknowledgements

We thank Karen Bustillo for helpful discussions. CO acknowledges support of a US Department of Energy Early Career Research Award. SEZ was supported by the National Science Foundation under STROBE Grant no. DMR 1548924. AB, BHS, and $\pyFDSTEM$ development are supported by the Toyota Research Institute. AR is supported by the 4D Data Distillery project, funded by the US Department of Energy. Work at the Molecular Foundry was supported by the Office of Science, Office of Basic Energy Sciences, of the U.S. Department of Energy under Contract No. DE-AC02-05CH11231. This research used resources of the National Energy Research Scientific Computing Center (NERSC), a U.S. Department of Energy Office of Science User Facility located at Lawrence Berkeley National Laboratory, operated under Contract No. DE-AC02-05CH11231.
