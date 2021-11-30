# Hessians and vibrational analysis

## Second energy derivatives in Hartree--Fock

In order to check whether a stationary point, i.e., a point on the potential energy surface with a vanishing gradient,
is a local minimum or not, the second derivatives of the energy with respect to nuclear displacement need to be calculated.
This matrix of all possible second derivatives is usually referred to as the **Hessian matrix**.
The Hessian can successively be diagonalized to obtain harmonic force constants, vibrational frequencies, and normal modes.
Second derivatives of the energy can, of course, be calculated numerically (based on either the numerical or analytical gradient,
which has the same advantages and disadvantages as in the case of the gradient itself), or analytically.
In the following, we will describe how the analytical Hessian is calculated at the level of Hartree--Fock (HF) theory.


The analytic gradient of the HF energy $E_{\text{HF}}$ with respect to a nuclear coordinate $\xi$, which was given in MO basis in the previous section,
can be equivalently written in AO basis as {cite}`Pople1979`
```{math}
:label: eq:HF_gradient_in_ao
  \frac{\mathrm{d} E_{\text{HF}}}{\mathrm{d} \xi} = \sum_{\mu \nu} P_{\mu \nu} h_{\mu \nu}^{\xi}
		+ \frac12 \sum_{\mu \nu \kappa \lambda} P_{\mu \nu} P_{\kappa \lambda} \langle \mu \kappa || \nu \lambda \rangle^\xi
		+ \sum_{\mu \nu} \omega_{\mu \nu} S_{\mu \nu}^{\xi} + \frac{\mathrm{d} V_{nn}}{\mathrm{d} \xi} \, ,
```
where the HF density matrix $\mathbf{P}$ and energy-weighted density matrix $\boldsymbol{\omega}$ in a real AO basis are defined as
```{math}
  P_{\mu \nu} &= \sum_{i} C_{\mu i} C_{\nu i} \\
  \omega_{\mu \nu} &= - \sum_{i} \varepsilon_i C_{\mu i} C_{\nu i}
```

Straightforward differentiation of the above equation {eq}`eq:HF_gradient_in_ao` with respect to another nuclear coordinate $\chi$ gives second derivatives of the HF energy:
```{math}
:label: eq:HF_Hessian_Pople
  \frac{\mathrm{d}^2 E_{\text{HF}}}{\mathrm{d} \chi \mathrm{d} \xi} &= \sum_{\mu \nu} P_{\mu \nu} h_{\mu \nu}^{\xi \chi} + \frac12 \sum_{\mu \nu \lambda \sigma} P_{\mu \nu} P_{\lambda \sigma} \langle \mu \lambda || \nu \sigma \rangle^{\xi \chi} + \sum_{\mu \nu} \omega_{\mu \nu} S_{\mu \nu}^{\xi \chi} + \frac{\mathrm{d}^2 V_{nn}}{\mathrm{d} \chi \mathrm{d} \xi} \\
  &+ \sum_{\mu \nu} \frac{\mathrm{d} P_{\mu \nu}}{\mathrm{d} \chi} h_{\mu \nu}^{\xi} + \sum_{\mu \nu \lambda \sigma} \frac{\mathrm{d} P_{\mu \nu}}{\mathrm{d} \chi} P_{\lambda \sigma} \langle \mu \lambda || \nu \sigma \rangle^{\xi} + \sum_{\mu \nu} \frac{\mathrm{d} \omega_{\mu \nu}}{\mathrm{d} \chi} S_{\mu \nu}^{\xi}
```

Now, calculation of the perturbed MO coefficients $C_{\mu p}^{\chi}$ can no longer be avoided {cite}`Pople1979`.

## Coupled-Perturbed Hartree--Fock

The perturbed MO coefficients $C_{\mu p}^{\chi}$ are expanded in the basis of the unperturbed ones {cite}`Pople1979`,
```{math}
:label: eq:perturbed_mo_coefficients
  C_{\mu p}^{\chi} &= \sum_{q} C_{\mu q} U_{q p}^{\chi} \, , \\
```
where the matrix $\mathbf{U}^\chi$ contains the unknown expansion coefficients. It can be shown that the occupied-occupied block of $\mathbf{U}^\chi$
is proportional to the partial derivative of the overlap matrix,
\begin{equation}
  U_{ij}^{\chi} = - \frac12 S_{ij}^{(\chi)} \, ,
\end{equation}
the virtual-virtual block is zero, $U_{ab}^\chi = 0$, and only the occupied-virtual block $U_{ai}^\chi$ needs to be determined.

The $U_{ai}^{\chi}$ are the solution of the **coupled-perturbed Hartree--Fock** (CPHF) equations, which can be written as {cite}`Deglmann2002`
\begin{align}
  (\varepsilon_i - \varepsilon_a) U_{ai}^{\chi} - 2 G_{ai}[U_{bj}^{\chi}] = R_{ai}^{\chi} \, ,
\end{align}
where the right-hand side (RHS), which depends on the perturbation $\chi$, is given by
\begin{equation}
  R_{ai}^{\chi} = F_{ai}^{(\chi)} - \varepsilon_{i} S_{ai}^{(\chi)} - G_{ai}[S_{jk}^{(\chi)}] \, ,
\end{equation}
and the matrices $\mathbf{G}[M_{rs}^{\chi}]$, for arbitrary $M_{rs}^{\chi}$, are defined as {cite}`Deglmann2002`
\begin{equation}
  G_{pq}[M_{rs}^{\chi}] = \sum_{\mu \nu \kappa \lambda} C_{\mu p} C_{\nu q} \langle \mu \kappa || \nu \lambda \rangle \sum_{rs} C_{\kappa r} C_{\lambda s} M_{rs}^{\chi} \, .
\end{equation}

It should be noted that in contrast to the first derivative of the MP2 or TDHF/TDA energy, where only a single orbital-response equation had to be solved,
a response equation needs to be solved for every perturbation $\chi$, i.e., the $3N$ Cartesian nuclear coordinates of a molecule containing $N$ atoms,
to obtain the second derivative of the HF energy.

Having solved for the CPHF coefficients $\mathbf{U}^{\chi}$, the perturbed density matrix can be calculated by using Eq. {eq}`eq:perturbed_mo_coefficients`
as
```{math}
:label: eq:perturbed_density
  \frac{\mathrm{d} P_{\mu \nu}}{\mathrm{d} \chi} &= \sum_{i} \frac{\mathrm{d}}{\mathrm{d} \chi} (C_{\mu i} C_{\nu i}) = \sum_{i} ( C_{\mu i}^\chi C_{\nu i} + C_{\mu i} C_{\nu i}^\chi ) \\
  &= - \frac12 \sum_{ij} (  C_{\mu i} C_{\nu j} S_{ij}^{(\chi)} + C_{\mu j} C_{\nu i} S_{ij}^{(\chi)} ) + \sum_{ia} ( C_{\mu a} C_{\nu i} + C_{\mu i} C_{\nu a} ) U_{ai}^{\chi} \\
```


As an alternative to Eq. {eq}`eq:HF_Hessian_Pople`, the second derivatives of the HF energy can then be written
in terms of the RHS $R_{ai}^{\chi}$, the CPHF coefficients $U_{ai}^{\xi}$, and partial derivatives of the Fock and overlap matrices as {cite}`Deglmann2002`
```{math}
:label: eq:HF_Hessian_Furche
  \frac{\mathrm{d}^2 E_{\text{HF}}}{\mathrm{d} \chi \mathrm{d} \xi} &= \sum_{\mu \nu} P_{\mu \nu} h_{\mu \nu}^{\chi \xi} + \frac12 \sum_{\mu \nu \lambda \sigma} P_{\mu \nu} P_{\lambda \sigma} \langle \mu \lambda || \nu \sigma \rangle^{\xi \chi} + \sum_{\mu \nu} \omega_{\mu \nu} S_{\mu \nu}^{\chi \xi} + \frac{\mathrm{d}^2 V_{nn}}{\mathrm{d} \chi \mathrm{d} \xi} \\
  &+ 2 \sum_{ia} R_{ai}^{\chi} U_{ai}^{\xi} - \sum_{ij} \Big( F_{ij}^{(\chi)} S_{ij}^{(\xi)} + F_{ij}^{(\xi)} S_{ij}^{(\chi)} - 2 \varepsilon_i S_{ij}^{(\chi)} S_{ij}^{(\xi)} - 2 G_{ij}[S_{kl}^{(\chi)}] S_{ij}^{(\xi)} \Big) 
```


%## Vibrational frequencies and normal modes
## Vibrational analysis

The following steps are carried out mostly by the geomeTRIC module {cite}`Wang2016`,
hence they are described in less detail.
For more details on the topic, the reader is referred to _Molecular Vibrations_
by Wilson, Decius and Cross {cite}`Wilson1980`.


### Cartesian and mass-weighted Hessian

The starting point for the vibrational analysis of molecules
is the Hessian matrix in Cartesian coordinates, $\mathbf{H}^{\text{Cart}}$,
calculated either numerically or analytically as described above.
Generally, the elements of $\mathbf{H}^{\text{Cart}}$ are given by second derivatives of the energy $E$
with respect to nuclear displacement,
```{math}
:label: eq:Hessian_matrix_elements
  H_{ij}^{\text{Cart}} = \bigg( \frac{\mathrm{d}^2 E}{\mathrm{d} \xi_i \mathrm{d} \xi_j} \bigg)_0 \, .
```
Hence, $\mathbf{H}^{\text{Cart}}$ is a $3N \times 3N$ matrix (where $N$ is the number of atoms),
and $\xi_1, \xi_2, \ldots, \xi_{3N}$ is used for the displacements
of the Cartesian coordiates, $\Delta x_1, \Delta y_1, \Delta z_1, \ldots, \Delta z_N$.
The '$0$' subscript of the parentheses refers to the equilibrium geometry of the atoms at which the derivatives are taken,
and that the first derivatives vanish.

As a first step, the Hessian is converted to _mass-weighted_ Cartesian coordinates (MWC),
$q_1 = \sqrt{m_1} \xi_1 = \sqrt{m_1} \Delta x_1$, $q_2 = \sqrt{m_1} \xi_2 = \sqrt{m_1} \Delta y_1$,
$\ldots$, $q_{3N} = \sqrt{m_N} \xi_{3N} = \sqrt{m_N} \Delta z_{N}$, where $m_i$ is the mass of atom $i$,
such that $\mathbf{H}^{\text{MWC}}$ is given by
```{math}
:label: eq:Hessian_mwc
  H_{ij}^{\text{MWC}} = \frac{H_{ij}^{\text{Cart}}}{\sqrt{m_i m_j}} = \bigg( \frac{\mathrm{d}^2 E}{\mathrm{d} q_i \mathrm{d} q_j} \bigg)_0 \, .
```
Diagonalizing this Hessian gives $3N$ eigenvalues are the fundamental frequencies of the molecule,
which still include the translation and rotational modes. However, these should be close to zero.


### Translating and rotating frame

In order to remove translational and rotational degrees of freedom,
one first determines the center of mass (COM) $\mathbf{R}^{\text{COM}}$ in the usual way,
```{math}
:label: eq:center_of_mass
  \mathbf{R}^{\text{COM}} = \frac{\sum_{K} m_{K} \mathbf{R}_K}{\sum_{K} \mathbf{R}_K} \, ,
```
where the sum runs over all atoms $K$, and the origin is then shifted to the COM,
$\mathbf{R}_{K}^{\text{COM}} = \mathbf{R}_K - \mathbf{R}^{\text{COM}}$.
Subsequently, one determines the inertia tensor and diagonalizes it
to obtain principal moments and axes of inertia.
Next, one needs to find the transformation
from mass-weighted Cartesian coordinates to a set of $3N$ coordinates,
where the molecule's translation and rotation are separated out,
leaving $3N - 6$ (or $3N-5$ for linear molecules) vibrational modes.

While the three vectors of length $3N$ corresponding to translation are simply given by $\sqrt{m_i}$
times the coordinate axis, the vectors corresponding to rotational motion of the atoms
are obtained from the coordinates of the atoms with respect to the COM
and the corresponding row of the matrix used to diagonalize the moment of inertia tensor.
This corresponds to internal coordinates in the Eckart frame {cite}`Eckart1934`.
In the next step, these vectors are normalized and a Gram--Schmidt orthogonalization
is carried out to create $N_\text{vib} = 3N-6$ (or $3N-5$) remaining vectors,
which are orthogonal to the five or six translational and rotational vectors.
Thus, one obtains a transformation matrix $\mathbf{D}$ which allows for the transformation
of the mass-weighted Cartesian coordinates $\mathbf{q}$ to internal coordinates
$\mathbf{S} = \mathbf{Dq}$, where translation and rotation have been projected out.
%Thus, one obtains a basis in which translation and rotation
%is projected out of the (mass-weighted) Cartesian coordinates.


### Hessian in internal coordinates and harmonic frequencies

Now the Hessian $\mathbf{H}^\text{MWC}$, which is still given in mass-weighted Cartesian coordinates,
is transformed the the internal coordinate system,
```{math}
:label: eq:Hessian_to_internal
  \mathbf{H}^{\text{Int}} = \mathbf{D}^\dagger \mathbf{H}^\text{MWC} \mathbf{D} \, ,
```
yielding a representation in $N_\text{vib}$ internal coordinates from the full $3N$ Cartesian coordinates.
The Hessian in internal coordinates $\mathbf{H}^{\text{Int}}$ is successively diagonalized,
```{math}
:label: eq:Hessian_internal_diagonalized
  \mathbf{L}^\dagger \mathbf{H}^{\text{Int}} \mathbf{L} = \mathbf{\Lambda} \, ,
```
where $\mathbf{\Lambda}$ is the diagonal matrix of $N_{\text{vib}}$ eigenvalues $\lambda_i$
which are related to the harmonic vibrational frequencies $\nu_i$
and $\mathbf{L}$ is the transformation matrix composed of the eigenvectors.

Finally, the eigenvalues $\lambda_i = 4 \pi^2 \nu_i^2$ can be converted from frequencies $\nu_i$
to wavenumbers $\tilde{\nu}_i$ in reciprocal centimeters by using the relationship
$\nu_i = c \tilde{\nu}_i$, where $c$ is the speed of light.
The wavenumbers are thus obtained from
```{math}
:label: eq:vibrations_wavenumbers
  \tilde{\nu}_i = \sqrt{\frac{\lambda_i}{4\pi^2 c^2}} \, ,
```
and successively appropriate conversion factors are applied to obtain the wavenumbers in cm$^{-1}$.


### Cartesian displacements, reduced masses, and force constants

The Cartesian normal modes $\mathbf{l}^{\text{Cart}}$ are obtained by combining Eqs. {eq}`eq:Hessian_to_internal`
and {eq}`eq:Hessian_internal_diagonalized` together with a diagonal matrix $\mathbf{M}$
defined by $M_{ii} = \frac{1}{\sqrt{m_i}}$ to undo the mass-weighting,
$\mathbf{l}^{\text{Cart}} = \mathbf{M D L}$, with the individual elements of this matrix being given by
```{math}
  l_{ij}^{\text{Cart}} = \sum_{k=1}^{3N} \frac{D_{ik} L_{kj}}{\sqrt{m_i}} \, .
```
The (normalized) column vectors of $\mathbf{l}^{\text{Cart}}$ correspond to the normal-mode displacements in Cartesian coordinates,
which are used for the calculation of spectroscopic properties as described below.

From the Cartesian normal modes $\mathbf{l}^{\text{Cart}}$, the reduced mass $\mu_i$ of vibration $i$
can be calculated as
```{math}
:label: eq:reduced_masses
  \mu_i = \frac{1}{\sum_{k=1}^{3N} \big( l_{ki}^{\text{Cart}} \big)^2} \, ,
```
and from those the corresponding force constants $k_i$ are calculated as
```{math}
:label: eq:force_constants
  k_i = 4 \pi^2 \tilde{\nu}_i^2 \mu_i \, ,
```
since $\tilde{\nu}_i = \frac{1}{2 \pi} \sqrt{\frac{k_i}{\mu_i}}$.
The force constants are then converted from atomic units to milli-dyne per ångström.



### Infrared intensities

In order to calculate intensities in the infrared (IR) spectrum,
the nuclear derivative of the electric dipole moment $\boldsymbol{\mu} = (\mu_x, \mu_y, \mu_z)$ is needed,
where each component $\mu$ can be decomposed into an electronic and a nuclear contribution, $\mu = \mu_{\text{e}} + \mu_{\text{n}}$.
The nuclear part of the dipole moment $\mu_{\text{n}}$ is simply given by the classical expression
```{math}
  \mu_{\text{n}} = \sum_{K} Z_K R_K
```
with the charge $Z_K$ and Cartesian coordinate $R_K$ of nucleus $K$.
The electronic part $\mu_{\text{e}}$ is calculed quantum-mechanically from the dipole moment integrals in AO basis,
$\mu_{\kappa \lambda} = \langle \phi_\kappa | q r | \phi_\lambda \rangle$,
where $q$ is the electron's charge and $r$ one of its coordinates, $r \in \{ x,y,z \}$,
and the one-particle density matrix $\mathbf{P}$,
```{math}
:label: eq:electric_dipole_moment
  \mu_{\text{e}} = \sum_{\kappa \lambda} P_{\lambda \kappa} \mu_{\kappa \lambda} \, .
```
The nuclear gradient of the dipole moment can again be calculated either numerically or analytically.
While the analytic derivative of the nuclear contribution $\mu_{\text{n}}$ is trivial,
the derivative of the electronic part {eq}`eq:electric_dipole_moment` is given by
```{math}
:label: eq:electronic_dipole_derivative
  \frac{\mathrm{d} \mu_{\text{e}}}{\mathrm{d} \chi} = \sum_{\kappa \lambda} \frac{\mathrm{d} P_{\lambda \kappa}}{\mathrm{d} \chi} \mu_{\kappa \lambda}
	+ \sum_{\kappa \lambda} P_{\lambda \kappa} \frac{\mathrm{d} \mu_{\kappa \lambda}}{\mathrm{d} \chi} \, ,
```
where the perturbed density from Eq. {eq}`eq:perturbed_density` and the derivatives of the dipole integrals are needed.

The IR transition dipole moment is then calculated by taking the dot product of the dipole moment gradient {eq}`eq:electronic_dipole_derivative`
with the Cartesian normal modes $\mathbf{l}^{\text{Cart}}$, and the IR intensity as the square norm of corresponding transition moment.
The intensities are successively converted to the unit of km mol$^{-1}$.
Raman intensities are calculated in an analogous manner, except that the nuclear derivative
of the electric-dipole polarizability $\boldsymbol{\alpha}$ is needed.





