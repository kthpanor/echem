(hessians:label)=
# Hessians

## Second Energy Derivatives in Hartree--Fock

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


## Vibrational Frequencies and Normal Modes

_write about mass-weighing the Hessian, diagonalizing it, obtaining harmonic frequencies
and Cartesian normal modes from it..._


### Infrared Intensities

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

While the derivative of the nuclear contribution $\mu_{\text{n}}$ is trivial,
the derivative of the electronic part {eq}`eq:electric_dipole_moment` is given by
```{math}
:label: eq:electronic_dipole_derivative
  \frac{\mathrm{d} \mu_{\text{e}}}{\mathrm{d} \chi} = \sum_{\kappa \lambda} \frac{\mathrm{d} P_{\lambda \kappa}}{\mathrm{d} \chi} \mu_{\kappa \lambda}
	+ \sum_{\kappa \lambda} P_{\lambda \kappa} \frac{\mathrm{d} \mu_{\kappa \lambda}}{\mathrm{d} \chi} \, ,
```
where the perturbed density from Eq. {eq}`eq:perturbed_density` and the derivatives of the dipole integrals are needed.






%If we work in AO basis:
%```{math}
%:label: eq:Lagrangian_in_ao
%\frac{\partial L}{\partial \xi}&=\frac{\partial E_\mathrm{HF}}{\partial \xi}+\sum_{p,q}\omega_{pq}\frac{\partial S_{pq}}{\partial \xi}\\
%&=\frac{\partial E_\mathrm{HF}}{\partial \xi}+\sum_{\mu\nu}\sum_{p,q}\omega_{pq}C_{\mu p}S^\xi_{\mu\nu}C_{\nu q}\\
%&=\frac{\partial E_\mathrm{HF}}{\partial \xi}+\sum_{\mu\nu}\omega_{\mu\nu}S^\xi_{\mu\nu}\\
%```
%and take the derivative:
%```{math}
%:label: eq:second_order_energy_derivative
%\frac{\mathrm{d}^2E}{\mathrm{d}\chi\mathrm{d}\xi}&=\frac{\mathrm{d}}{\mathrm{d}\chi}\frac{\partial L}{\partial \xi}=\frac{\mathrm{d}}{\mathrm{d}\chi}\left(\frac{\partial E_\mathrm{HF}}{\partial \xi}+\sum_{\mu\nu}\omega_{\mu\nu}S^\xi_{\mu\nu}\right)\\
%&=\frac{\mathrm{d}}{\mathrm{d}\chi}\frac{\partial E_\mathrm{HF}}{\partial \xi} + \sum_{\mu\nu}\frac{\mathrm{d}}{\mathrm{d}\chi}\left(\omega_{\mu\nu}S^\xi_{\mu\nu}\right)\\
%&=\frac{\mathrm{d}}{\mathrm{d}\chi}\frac{\partial E_\mathrm{HF}}{\partial \xi} + \sum_{\mu\nu}\left(\frac{\mathrm{d}\omega_{\mu\nu}}{\mathrm{d}\chi}S^\xi_{\mu\nu}+\omega_{\mu\nu}\frac{\mathrm{d}}{\mathrm{d}\chi}S^\xi_{\mu\nu}\right)\\
%&=\frac{\mathrm{d}}{\mathrm{d}\chi}\frac{\partial E_\mathrm{HF}}{\partial \xi} + \sum_{\mu\nu}\left(\frac{\mathrm{d}\omega_{\mu\nu}}{\mathrm{d}\chi}S^\xi_{\mu\nu}+\omega_{\mu\nu}S^{\xi\chi}_{\mu\nu}\right)\\
%```
%By plugging in the expression for the partial derivative of the HF energy in AO basis, we get:
%```{math}
%:label: eq:second_order_explicit
%\frac{\mathrm{d}^2E}{\mathrm{d}\chi\mathrm{d}\xi}&=\frac{\mathrm{d}}{\mathrm{d}\chi}\left( \sum_{\mu \nu} P_{\mu \nu} h_{\mu \nu}^{\xi} + \frac12 \sum_{\mu \nu \lambda \sigma} P_{\mu \nu} P_{\lambda \sigma} \langle \mu \lambda || \nu \sigma \rangle^{\xi} + \frac{\partial V_{nn}}{\partial \xi}\right) + \sum_{\mu\nu}\left(\frac{\mathrm{d}\omega_{\mu\nu}}{\mathrm{d}\chi}S^\xi_{\mu\nu}+\omega_{\mu\nu}S^{\xi\chi}_{\mu\nu}\right)\\
%&= \sum_{\mu \nu} \left( \frac{\mathrm{d}P_{\mu \nu}}{\mathrm{d}\chi}  h_{\mu \nu}^{\xi} + P_{\mu \nu}h_{\mu\nu}^{\xi\chi}  \right) + \frac12 \sum_{\mu \nu \lambda \sigma} \left( \frac{\mathrm{d}P_{\mu \nu}}{\mathrm{d}\chi}  P_{\lambda \sigma} \langle \mu \lambda || \nu \sigma \rangle^{\xi} + P_{\mu \nu}\frac{\mathrm{d}P_{\lambda \sigma}}{\mathrm{d}\chi}   \langle \mu \lambda || \nu \sigma \rangle^{\xi}+P_{\mu \nu}P_{\lambda \sigma} \langle \mu \lambda || \nu \sigma \rangle^{\xi\chi} \right) \\
%&+ \frac{\partial^2 V_{nn}}{\partial \xi\partial \chi} + \sum_{\mu\nu}\left(\frac{\mathrm{d}\omega_{\mu\nu}}{\mathrm{d}\chi}S^\xi_{\mu\nu}+\omega_{\mu\nu}S^{\xi\chi}_{\mu\nu}\right)\\
%&= \sum_{\mu \nu} \left( \frac{\mathrm{d}P_{\mu \nu}}{\mathrm{d}\chi}  h_{\mu \nu}^{\xi} + P_{\mu \nu}h_{\mu\nu}^{\xi\chi}  \right) + \frac12 \sum_{\mu\nu\lambda\sigma}P_{\mu \nu}P_{\lambda \sigma} \langle \mu \lambda || \nu \sigma \rangle^{\xi\chi} + \sum_{\mu \nu \lambda \sigma} \frac{\mathrm{d}P_{\mu \nu}}{\mathrm{d}\chi}  P_{\lambda \sigma} \langle \mu \lambda || \nu \sigma \rangle^{\xi} +\frac{\partial^2 V_{nn}}{\partial \xi\partial \chi}\\
%&+ \sum_{\mu\nu}\left(\frac{\mathrm{d}\omega_{\mu\nu}}{\mathrm{d}\chi}S^\xi_{\mu\nu}+\omega_{\mu\nu}S^{\xi\chi}_{\mu\nu}\right)
%```
