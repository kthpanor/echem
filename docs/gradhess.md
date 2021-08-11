# Gradients and Hessians
To add: general aspects about gradients and Hessians, numerical vs. analytical, related properties.

## Analytical Gradients
To add: Introduction, approaches to derive analytical gradients; Lagrangian approach.

To determine the energy gradient, we first need to identify the non-variational components of the energy functional. In the case of DFT (and HF), these are the MO coefficients ($\mathbf{C}$), which exhibit an implicit dependence on the nuclear coordinates when atom-centred basis functions such as Gaussian or Slater-type orbitals are used. Considering this implicit dependence, the total energy derivative with respect to a particular nuclear coordinate $\xi$ is obtained via the chain rule {cite}`Rehn2015,Levchenko2005`:
%
```{math}
:label: eq:energy_functional_general
\frac{dE}{d\xi}=\frac{\partial E}{\partial \xi}+\frac{\partial E}{\partial\mathbf{C}}\frac{d\mathbf{C}}{d\xi} \, .
```
%
Here, $\mathbf{C}$ is the molecular orbital matrix which transforms from a set of atomic orbitals $\{\phi_\mu\}$ to a set of molecular orbitals $\{\psi_p\}$ via $\psi_p=\sum_\mu C_{\mu p}\phi_\mu$. The first term, $\partial E/\partial\xi$, is the Hellman--Feynman contribution which describes the explicit dependence of the energy on the nuclear coordinate $\xi$ through the nuclear-electron and nuclear-nuclear interaction terms of the Hamiltonian {cite}`Levchenko2005_thesis, Helgaker1988_analytical`. The second term stems from the implicit dependence of the energy on $\xi$ due to the fact that the molecular orbitals are expanded in a finite atomic-centred basis set {cite}`Helgaker1988_analytical` 
%This term would vanish if a basis set independent on atomic positions, for example a plane wave basis, were used.

It may seem at first surprising that the derivative ${\partial E}/{\partial\mathbf{C}}$ has to be computed. If the MO coefficients are obtained variationally for a specified molecular geometry, how is it that this derivative is not zero? The key to this conundrum lies in the phrase "for a specified molecular geometry". Since the DFT energy and density are constructed using a constrained LCAO parametrization, if we perform a nuclear displacement, the "old" MO coefficients no longer correspond to the minimum energy and must be re-optimized. Thus the partial derivative with respect to the MO coefficients, as well as the derivative of the MO coefficients with respect to $\xi$ are required. The explicit computation of the latter is complicated, but can be avoided by making use of a new functional, the Lagrangian, for which the partial derivative $\partial L / \partial \mathbf{C}$ is zero and constrained to the DFT/HF configuration space by construction {cite}`Levchenko2005, Helgaker2014`:
%ensuring that we re still in the H  Thus, the derivative
%Thus, even though the local points of the PES are variational with respect to the MO coefficients, the "global" energy functional is not. We must therefore, include the partial derivative with respect to the MO coefficients, as well as the derivative of the MO coefficients with respect to $\xi$, which is quite complicated to compute. However, we can avoid explicitly computing $\mathrm{d}\mathbf{C}/\mathrm{d}\xi$ by using a trick. The idea is to create a new functional, the Lagrangian ($L$), which is by construction equivalent to the energy functional, but for which the partial derivative $\partial L / \partial \mathbf{C}$ is zero \cite{Levchenko2005, Helgaker2014}:
%
```{math}
:label: eq:Lagrangian_general
L(\mathbf{C},\boldsymbol{\Lambda})=E(\mathbf{C})+\boldsymbol{\Lambda}f_c(\mathbf{C}) \, ,
```
%
where $\boldsymbol{\Lambda}$ are a set of undetermined Lagrange multipliers and $f_c(\mathbf{C})=0$ define the constraints for the non-variational parameters $\mathbf{C}$. These constraints ensure that we are moving only in the DFT/HF configuration space, rather than in the infinite space of all orthogonally equivalent combinations of orbital bases {cite}`Helgaker1988_analytical`.
%From Helgaker and Jorgensen:"At each geometry we have a set of AOs  from  which an infinite set of orthogonally equivalent orbital bases can be constructed. As the geometry changes we must pick out exactly one of these orbital bases at each geometry X. In this way an orthogonal  orbital connection is established. (A connection is called orthogonal if it  preserves orthonormality between the orbitals.) We further  require that the connection is continuous and differentiable. One may also wish to impose an additional requirement on the connection, namely that it  is translationally and rotationally  invariant. This may seem to be a trivial requirement. However, a connection is conveniently defined  in terms of atomic Cartesian displacements rather  than in terms of a set of nonredundant internal  coordinates. This implies that each  molecular geometry  may be described  in an infinite number of translationally and rotationally  equivalent ways"

By using the Lagrangian, we have shifted the difficult problem of computing ${\mathrm{d} \mathbf{C}}/{\mathrm{d}\xi}$ to the much simpler problem of determining the unknown Lagrange multipliers which satisfy $\partial L / \partial \mathbf{C}=0$. Equations for these are derived by imposing that the explicit form of $\partial L / \partial \mathbf{C}$ is zero {cite}`Helgaker2014`. Once the Lagrange multipliers have been obtained, the total derivative of the energy functional with respect to the nuclear coordinate $\xi$ can be computed as:
\begin{equation}
\frac{dE}{d\xi}=\frac{\partial L}{\partial\xi}\, .
\end{equation}

### Ground state
#### Hartree--Fock
As an example, we will derive in the following the analytical expression for the Hartree-Fock energy gradient. The DFT gradient can be derived in a similar way, replacing the exact exchange integrals with the corresponding exchange--correlation functional contributions. The exchange-correlation functional contribution to the DFT energy and its molecular gradient is evaluated via numerical integration. Thus, the molecular gradient includes grid point weight contributions, which arise from the explicit dependence of the grid partitioning function on the molecular geometry. Neglecting these contributions to the molecular gradient leads to the breakdown of rotation--translation invariance of the molecular gradient. Despite this, if a fine integration grid is used in practical calculations, grid point weight contribution to the molecular gradient can be safely neglected.  

We are interested to compute the derivative of the electronic energy of the ground state $\ket{0}$, described at the HF level of theory, with respect to a nuclear coordinate $\xi$. 

 ```{math}
:label: eq:HF_energy_fct
E_\mathrm{HF}=\sum_{i}f_{ii}-\frac{1}{2}\sum_{ij}\braket{ij||ij}\, ,
```
where $f_{ii}$ are Fock matrix elements, $\braket{ij||ij}$ are anti-symmetrized two electron integrals in Physicist's notation, and indices $i,j$ denote occupied molecular orbitals.

The Lagrangian for this energy functional is constructed as follows:

 ```{math}
:label: eq:HF_Lagrangian
L(\mathbf{C},\boldsymbol{\Lambda},\boldsymbol{\Omega})=E_\mathrm{HF}+\sum_{p,q}\lambda_{pq}\left(f_{pq}-\delta_{pq}\epsilon_p\right)+\sum_{p,q}\omega_{pq}\left(S_{pq}-\delta_{pq}\right) \, ,
```
with $\boldsymbol{\Lambda}=\{\lambda_{pq}\}$ and $\boldsymbol{\Omega}=\{\omega_{pq}\}$ the Lagrange multipliers, $\mathbf{F}=\{f_{pq}\}$ the Fock matrix, $\{\epsilon_p\}$ the Hartree-Fock orbital energies, and $\mathbf{S}=\{S_{pq}\}$ the overlap matrix. Here, the conditions for the Lagrange multipliers $\{\lambda_{pq}\}$ and $\{\omega_{pq}\}$ ensure that the Fock matrix is diagonal and, respectively, the overlap matrix is unity for the HF state.

To calculate the total derivative of the energy with respect to $\xi$ we must now only calculate the partial derivative of the Lagrangian with respect to the same variable:
```{math}
:label: eq:partial_L
\frac{\partial L}{\partial \xi}=\frac{\partial E_\mathrm{HF}}{\partial \xi}+\sum_{p,q}\lambda_{pq}\frac{\partial f_{pq}}{\partial \xi}+\sum_{p,q}\omega_{pq}\frac{\partial S_{pq}}{\partial \xi}\, .
```
We are interested to re-write the above derivative of the Lagrangian in terms of effective one- and two-particle density matrices. For this purpose, we express the energy in terms of the one- $\{\gamma_{pq}\}$ and two-particle $\{\Gamma_{pqrs}\}$ density matrices:
```{math}
:label: eq:EHF_DM
E_\mathrm{HF}=\sum_{p,q}f_{pq}\gamma_{pq}+\frac{1}{4}\sum_{pqrs}\Gamma_{pqrs}\braket{pq||rs}\, .
```

With this definition, Eq. {eq}`eq:partial_L` becomes:
```{math}
:label: eq:partial_L_final
\frac{\partial L}{\partial \xi}&=\sum_{p,q}(\lambda_{pq}+\gamma_{pq}) f^\xi_{pq}+\frac{1}{4}\sum_{pqrs}\Gamma_{pqrs}\braket{pq||rs}^\xi+\sum_{p,q}\omega_{pq}S^\xi_{pq}\nonumber\\
&=\sum_{p,q}(\lambda_{pq}+\gamma_{pq}) h^\xi_{pq}+\sum_{p,q}(\lambda_{pq}+\gamma_{pq})\sum_{i}\braket{pi||qi}^\xi+\frac{1}{4}\sum_{pqrs}\Gamma_{pqrs}\braket{pq||rs}^\xi+\sum_{p,q}\omega_{pq}S^\xi_{pq}\,,
```
where $h_{pq}$ represents a matrix element of the core-Hamiltonian operator. The superscript $\xi$ indicates a partial derivative with respect to variable $\xi$, explicitly {cite}`Levchenko2005`:
```{math}
:label: eq:hpq
h_{pq}^{\xi}&=\frac{\partial h_{pq}}{\partial \xi}=\sum_{\mu,\nu}C_{\mu p}h^{\xi}_{\mu\nu}C_{\nu q}\, ,\\
h^{\xi}_{\mu\nu}&=\bra{\phi_\mu}\frac{\partial \hat{h}}{\partial \xi}\ket{\phi_\nu}+\bra{\frac{\partial\phi_\mu}{\partial \xi}}\hat{h}\ket{\phi_\nu}+\bra{\phi_\mu}\hat{h}\ket{\frac{\partial\phi_\nu}{\partial\xi}}\, ,\\
```
```{math}
:label: eq:braket_pqrs
\braket{pq||rs}^{\xi}&=\frac{\partial \braket{pq||rs}}{\partial\xi}=\sum_{\mu,\nu,\theta,\sigma}C_{\mu p}C_{\nu q}\braket{\phi_{\mu}\phi_{\nu}||\phi_{\theta}\phi_{\sigma}}^{\xi}C_{\theta r}C_{\sigma s}\, ,\\
\braket{\phi_{\mu}\phi_{\nu}||\phi_{\theta}\phi_{\sigma}}^{\xi}&=\braket{\frac{\partial\phi_{\mu}}{\partial\xi}\phi_{\nu}||\phi_{\theta}\phi_{\sigma}}+\braket{\phi_{\mu}\frac{\partial\phi_{\nu}}{\partial\xi}||\phi_{\theta}\phi_{\sigma}}+\braket{\phi_{\mu}\phi_{\nu}||\frac{\partial\phi_{\theta}}{\partial\xi}\phi_{\sigma}}\nonumber\\
&+\braket{\phi_{\mu}\phi_{\nu}||\phi_{\theta}\frac{\partial\phi_{\sigma}}{\partial \xi}}\,\\
```
```{math}
:label: eq:Spq
S_{pq}^\xi &=\frac{\partial S_{pq}}{\partial\xi}=\sum_{\mu,\nu}C_{\mu p}S_{\mu\nu}^\xi C_{\nu q}\, ,\\
S_{\mu\nu}^\xi &=\braket{\frac{\partial\phi_\mu}{\partial\xi}|\phi_\nu}+\braket{\phi_\mu|\frac{\partial\phi_\nu}{\partial\xi}} \,.
```

We also made use of the definition of the Fock matrix {cite}`Szabo2012`: 
```{math}
:label: eq:FockMatEl
f_{pq}=h_{pq}+\sum_i\braket{pi||qi}\, .
```
Equation {eq}`eq:partial_L_final` is the working equation to compute the partial derivative of the Lagrangian with respect to variable $\xi$, and implicitly the equation for the total derivative of the energy with respect to the same variable. Two ingredients are required to calculate $\partial L/\partial\xi$: (1) the derivatives of the core-Hamiltonian matrix, of the anti-symmetrized two-electron integrals, and of the overlap matrix, and (2) finding the density matrices ($\gamma_{pq}$, $\Gamma_{pqrs}$) and Lagrange multipliers ($\lambda_{pq}$, $\omega_{pq}$).

By comparing Eq. {eq}`eq:EHF_DM` to Eq. {eq}`eq:HF_energy_fct`, we can immediately identify the density matrices:
```{math}
:label: eq:HF_gamma_ij
\gamma_{ij}=\delta_{ij} \, , 
```
```{math}
:label: eq:HF_Gamma_ijkl
\Gamma_{ijkl}=-2\delta_{ik}\delta_{jl}\, .
```
Equations for the $\{\lambda_{pq}\}$ and $\{\omega_{pq}\}$ multipliers are obtained by imposing the Lagrangian to be stationary with respect to the orbital transformation matrix $\{C_{\mu t}\}$:
```{math}
:label: eq:OrbitalRspCond
\frac{\partial L}{\partial C_{\mu t}}=0\, ,
```
or the programmable version {cite}`Levchenko2005,Rehn2019` :
```{math}
:label: eq:OrbitalRspCond_program
\sum_{\mu}C_{\mu u}\frac{\partial L}{\partial C_{\mu t}}=0\, ,
```
To calculate the partial derivative of the Lagrangian with respect to $C_{\mu t}$, we will need the following three expressions:
```{math}
:label: eq:derivFock
\sum_{\mu}C_{\mu u}\frac{\partial f_{pq}}{\partial C_{\mu t}}&=\sum_{\mu}C_{\mu u}\left(\frac{\partial h_{pq}}{\partial C_{\mu t}}+\sum_i\frac{\partial \braket{pi||qi}}{\partial C_{\mu t}}\right)\nonumber\\
&=h_{uq}\delta_{pt}+h_{pu}\delta_{qt}+\sum_i\braket{ui||qi}\delta_{pt}+\sum_i\braket{pi||ui}\delta_{qt}\nonumber\\
&+\sum_i\braket{pu||qi}\delta_{it}+\sum_i\braket{pi||qu}\delta_{it}\nonumber \\
&=f_{uq}\delta_{pt}+f_{pu}\delta_{qt}+\braket{pu||qt}\delta_{t\epsilon_o}+\braket{pt||qu}\delta_{t\epsilon_o} \, ,\label{eq:derivFock}
```
```{math}
:label: eq:derivTwoEl
\sum_{\mu}C_{\mu u}\frac{\partial \braket{pq||rs}}{\partial C_{\mu t}}=\braket{uq||rs}\delta_{pt}+\braket{pu||rs}\delta_{qt}+\braket{pq||us}\delta_{rt}+\braket{pq||ru}\delta_{st} \, ,
```
```{math}
:label: eq:derivOverlap
\sum_{\mu}C_{\mu u}\frac{\partial S_{pq}}{\partial C_{\mu t}}=S_{uq}\delta_{pt}+S_{pu}\delta_{qt}\, , \label{eq:derivOverlap}
```
where we have used the definitions of the Fock matrix, two-electron integrals and overlap matrix in terms of the orbital transformation matrix $\{C_{\mu p}\}$ -- see Eqs. {eq}`eq:hpq`, {eq}`eq:braket_pqrs`, and {eq}`eq:Spq`. The Kronecker delta $\delta_{t\epsilon_\mathrm{o}}$ is equal to one if $t$ is an occupied orbital and is zero otherwise.

Using the Lagrangian expressed in terms of the density matrices:
```{math}
:label: eq:L_DMs
L=\sum_{p,q}\gamma_{pq}f_{pq}+\sum_{p,q}\lambda_{pq}(f_{pq}-\delta_{pq}\epsilon_p)+\frac{1}{4}\sum_{pqrs}\Gamma_{pqrs}\braket{pq||rs}+\sum_{p,q}\omega_{pq}(S_{pq}-\delta_{pq})\,,
```

the partial derivative of the Lagrangian with respect to $C_{\mu p}$ can be written as:
```{math}
:label: eq:OrbRsp1
\sum_{\mu}C_{\mu u}\frac{\partial L}{\partial C_{\mu t}}&=\sum_{p,q}\left(\gamma_{pq}+\lambda_{pq}\right)\left[f_{uq}\delta_{pt}+f_{pu}\delta_{qt}+\braket{pu||qt}\delta_{t\epsilon_o}+\braket{pt||qu}\delta_{t\epsilon_o}\right]\nonumber \\
&+\frac{1}{4}\sum_{p,q,r,s}\Gamma_{pqrs}\left[\braket{uq||rs}\delta_{pt}+\braket{pu||rs}\delta_{qt}+\braket{pq||us}\delta_{rt}+\braket{pq||ru}\delta_{st}\right]\nonumber\\
&+\sum_{p,q}\omega_{pq}\left(S_{uq}\delta_{pt}+S_{pu}\delta_{qt}\right)\, .
```
By using the conditions $f_{pq}=\epsilon_p\delta_{pq}$ and $S_{pq}=\delta_{pq}$, Eq. {eq}`eq:OrbRsp1` becomes:
```{math}
:label: eq:OrbRsp2
\sum_{\mu}C_{\mu u}\frac{\partial L}{\partial C_{\mu t}}&=\sum_{p,q}\left(\gamma_{pq}+\lambda_{pq}\right)\left[\epsilon_u\delta_{uq}\delta_{pt}+\epsilon_u\delta_{pu}\delta_{qt}+\braket{pu||qt}\delta_{t\epsilon_o}+\braket{pt||qu}\delta_{t\epsilon_o}\right]\nonumber \\
&+\frac{1}{4}\sum_{p,q,r,s}\Gamma_{pqrs}\left[\braket{uq||rs}\delta_{pt}+\braket{pu||rs}\delta_{qt}+\braket{pq||us}\delta_{rt}+\braket{pq||ru}\delta_{st}\right]\nonumber\\
&+\sum_{p,q}\omega_{pq}\left(\delta_{uq}\delta_{pt}+\delta_{pu}\delta_{qt}\right) \nonumber \\
&=2\left(\gamma_{ut}+\lambda_{ut}\right)\epsilon_u+2\sum_{p,q}\left(\gamma_{pq}+\lambda_{pq}\right)\braket{pu||qt}\delta_{t\epsilon_o}+\sum_{p,q,r}\Gamma_{tpqr}\braket{up||qr}+2\omega_{ut}\, ,
```
where we have used $\gamma_{pq}=\gamma_{qp}$, $\braket{pu||qt}=\braket{qt||pu}$, $\Gamma_{pqrs}=\Gamma_{qpsr}=\Gamma_{srpq}$ (real orbitals), and we have imposed that $\lambda_{pq}=\lambda_{qp}$, and $\omega_{pq}=\omega_{qp}$ (symmetric representation). Some of the indices have been renamed.

To obtain equations for the orbital response Lagrange multipliers, we first have to decouple $\{\lambda\}$ from $\{\omega\}$ by taking the difference
```{math}
:label: eq:decouple
\sum_{\mu}C_{\mu u}\frac{\partial L}{\partial C_{\mu t}}-\sum_{\mu}C_{\mu t}\frac{\partial L}{\partial C_{\mu u}}&=2(\gamma_{ut}+\lambda_{ut})(\epsilon_u-\epsilon_t)\nonumber\\&+2\sum_{p,q}(\gamma_{pq}+\lambda_{pq})\left(\braket{pu||qt}\delta_{t\epsilon_o}-\braket{pt||qu}\delta_{u\epsilon_o}\right)\nonumber\\
&+\,\,\,\,\sum_{p,q,r}\left(\Gamma_{tpqr}\braket{up||qr}-\Gamma_{upqr}\braket{tp||qr}\right).
```
The system of equations for $\{\lambda\}$ is then obtained by choosing $u$ and $t$ from different orbital spaces in the following equation,
```{math}
:label: eq:OrbRspEq
&2(\gamma_{ut}+\lambda_{ut})(\epsilon_u-\epsilon_t)+2\sum_{p,q}(\gamma_{pq}+\lambda_{pq})\left(\braket{pu||qt}\delta_{t\epsilon_o}-\braket{pt||qu}\delta_{u\epsilon_o}\right)\nonumber\\
&+\sum_{p,q,r}\left(\Gamma_{tpqr}\braket{up||qr}-\Gamma_{upqr}\braket{tp||qr}\right)=0\,.
```
Once $\{\lambda\}$ is determined, $\{\omega\}$ can be calculated in a similar way, using the following equation
```{math}
:label: eq:omega
2\left(\gamma_{ut}+\lambda_{ut}\right)\epsilon_u+2\sum_{p,q}\left(\gamma_{pq}+\lambda_{pq}\right)\braket{pu||qt}\delta_{t\epsilon_o}+\sum_{p,q,r}\Gamma_{tpqr}\braket{up||qr}+2\omega_{ut}=0
```
If we explicitly write the equations for different blocks of $\{\lambda\}$, we find that they are all zero. This simplifies the equation for $\{\omega\}$ to:
```{math}
:label: eq:omega_HF
2\gamma_{ut}\epsilon_u+2\sum_{p,q}\gamma_{pq}\braket{pu||qt}\delta_{t\epsilon_o}+\sum_{p,q,r}\Gamma_{tpqr}\braket{up||qr}+2\omega_{ut}=0
```
The only non-zero block of $\{\omega\}$ is the occupied-occupied block:
```{math}
:label: eq:omega_hf_oo
&u=i,\, t=j \nonumber\\
&2\gamma_{ij}\,\epsilon_i+2\sum_{p,q}\gamma_{pq}\braket{pi||qj}\delta_{j\epsilon_o}+\sum_{p,q,r}\Gamma_{jpqr}\braket{ip||qr}+2\omega_{ij} = 0 \nonumber \\
&\Leftrightarrow
2\delta_{ij}\epsilon_i+2\sum_{k,l}\delta_{kl}\braket{ki||lj}+\sum_{klm}(-2\delta_{jl}\delta_{km})\braket{ik||lm}+2\omega_{ij} =0\nonumber \\
&\Leftrightarrow
\omega_{ij}=-\epsilon_i\delta_{ij},
```
Using the expressions for the density matrices and the non-zero Lagrange multipliers, we get the following expression for the electronic energy derivative:
```{math}
:label: eq:HF_grad_final
\frac{\mathrm{d}E}{\mathrm{d}\xi}&=\sum_{i,j}\delta_{ij} h^\xi_{ij}+\sum_{i,j}\delta_{ij}\sum_{k}\braket{ik||jk}^\xi-\frac{1}{2}\sum_{i,j,k,l}\delta_{ik}\delta_{jl}\braket{ij||kl}^\xi-\sum_{i,j}\epsilon_{i}\delta_{ij}S^\xi_{ij}\nonumber \\
&=\sum_{i}h^\xi_{ii}+\frac{1}{2}\sum_{i,j}\braket{ij||ij}^\xi-\sum_{i}\epsilon_{i}S^\xi_{ii}
```
Finally, the derivative of the total energy is obtained by adding the trivial contribution from the nuclear repulsion energy term {cite}`Szabo2012`.

#### DFT

#### MP2
In the case of M{\o}ller--Plesset theory, the energy functional has additional non-variational parameters that have to be considered when computing the gradient. These are the so-called t-amplitudes $\mathbf{T}$, so the corresponding term which has to be determined is called amplitude response. 
```{math}
:label: eq:energy_functional_MP
\frac{dE}{d\xi}=\frac{\partial E}{\partial \xi}+\frac{\partial E}{\partial\mathbf{C}}\frac{d\mathbf{C}}{d\xi}+\frac{\partial E}{\partial\mathbf{T}}\frac{d\mathbf{T}}{d\xi} \, .
```
The analytic expression for the MP energy gradient is obtained in a very similar way as we have done for the SCF ground state. The difference is that the Lagrangian contains additional Lagrange multipliers and constraints for the t-amplitudes. After obtaining the corresponding amplitude response Lagrange multipliers, these additional contributions will be written in terms of one- and two-particle density matrices, exactly as the total energy. Let's illustrate the procedure for the second order M{\o}ller--Plesset theory, MP2. At this level of theory, the total energy functional is written as:
```{math}
:label: eq:MP2_energy_fct
E_0=E_\mathrm{HF}+E_\mathrm{MP2}=\sum_{i}f_{ii}-\frac{1}{2}\sum_{ij}\braket{ij||ij}-\frac{1}{4}\sum_{i,j,a,b}\braket{ij||ab}t_{ijab}, \label{eq:MP2_energy_fct}
```
where
```{math}
:label: eq:MP2_tamplitudes
t_{ijab} = \frac{\braket{ij||ab}}{\epsilon_a+\epsilon_b-\epsilon_i-\epsilon_j} \label{eq:MP2_tamplitudes}
```
are the MP2 t-amplitudes.

The Lagrangian corresponding to this energy functional is:
```{math}
:label: eq:MP2_Lagrangian
L(\mathbf{C}, \mathbf{T},\boldsymbol{\Lambda},\boldsymbol{\Omega},\mathbf{\tilde{T}})=E_\mathrm{0}+\sum_{p,q}\lambda_{pq}\left(f_{pq}-\delta_{pq}\epsilon_p\right)+\sum_{p,q}\omega_{pq}\left(S_{pq}-\delta_{pq}\right)+\sum_{i,j,a,b}\tilde{t}_{ijab}f(t_{ijab})\,.
```
Here, $\tilde{T}=\{\tilde{t}_{ijab}\}$ are the amplitude response Lagrange multipliers and $f(t_{ijab})=0$ is the constraint. For MP2, this is:
```{math}
:label: eq:t_condition
f(t_{ijab})=t_{ijab}\left(\epsilon_a+\epsilon_b-\epsilon_i-\epsilon_j\right)-\braket{ij||ab}\,  .
```

The amplitude response Lagrange multipliers are determined by imposing the Lagrangian to be stationary with respect to the t-amplitudes.
```{math}
:label: eq:amplitude_rsp_general
\frac{\partial L}{\partial \mathbf{T}}=0\,.
```
Replacing $L$ and $\mathbf{T}$ in the equation above with the corresponding MP2 expressions, we get:
```{math}
:label: eq:t_cond_explicit
\frac{\partial L}{\partial t_{ijab}}&=\frac{\partial E_\mathrm{MP2}}{\partial t_{ijab}}+\sum_{k,l,c,d}\frac{\partial \tilde{t}_{klcd}f(t_{klcd})}{\partial t_{ijab}}\nonumber\\
&=-\frac{1}{4}\sum_{k,l,c,d}\braket{kl||cd}\frac{\partial t_{klcd}}{\partial t_{ijab}}+\sum_{k,l,c,d}\tilde{t}_{klcd}\frac{\partial t_{klcd}}{\partial t_{ijab}}\left(\epsilon_c+\epsilon_d-\epsilon_k-\epsilon_l\right)\nonumber\\
&=-\braket{ij||ab}+4\tilde{t}_{ijab}\left(\epsilon_a+\epsilon_b-\epsilon_i-\epsilon_j\right)\, , \label{eq:t_cond_explicit}
```
From Eqs. {eq}`eq:amplitude_rsp_general` and {eq}`eq:t_cond_explicit`, it follows that,
```{math}
:label: eq:ampl_rsp_multipliers
\tilde{t}_{ijab} = \frac{1}{4}\frac{\braket{ij||ab}}{\epsilon_a+\epsilon_b-\epsilon_i-\epsilon_j}=\frac{1}{4}t_{ijab}.
```
From here, we can follow the same procedure as we did for the SCF gradient. We first rewrite the Lagrangian in terms of one- and two-particle density matrices:
```{math}
:label: eq:L_DMs_MP
L&=\sum_{p,q}\gamma_{pq}f_{pq}+\sum_{p,q}\lambda_{pq}(f_{pq}-\delta_{pq}\epsilon_p)+\frac{1}{4}\sum_{pqrs}\Gamma_{pqrs}\braket{pq||rs}\nonumber\\
&+\sum_{p,q}\omega_{pq}(S_{pq}-\delta_{pq})+\sum_{i,j,a,b}\tilde{t}_{ijab}\left[t_{ijab}\left(\epsilon_a+\epsilon_b-\epsilon_i-\epsilon_j\right)-\braket{ij||ab}\right]\nonumber\\
&=\sum_{p,q}\gamma_{pq}f_{pq}+\sum_{p,q}\lambda_{pq}(f_{pq}-\delta_{pq}\epsilon_p)+\frac{1}{4}\sum_{pqrs}\Gamma_{pqrs}\braket{pq||rs}\nonumber\\
&+\sum_{p,q}\omega_{pq}(S_{pq}-\delta_{pq})+\sum_{p,q}\gamma^\mathrm{A}_{pq}f_{pq}+\frac{1}{4}\sum_{pqrs}\Gamma^\mathrm{A}_{pqrs}\braket{pq||rs}\,,
```
where we have written also the amplitude contribution in terms of one- and two-particle density matrices, $\gamma^\mathrm{A}_{pq}$ and $\Gamma^\mathrm{A}_{pqrs}$ respectively. Denoting $\gamma'=\gamma+\gamma^\mathrm{A}$ and $\Gamma'=\Gamma+\Gamma^\mathrm{A}$, the Lagrangian becomes:
```{math}
:label: eq:final_L_DMs_MP
L=\sum_{p,q}\gamma'_{pq}f_{pq}+\sum_{p,q}\lambda_{pq}(f_{pq}-\delta_{pq}\epsilon_p)+\frac{1}{4}\sum_{pqrs}\Gamma'_{pqrs}\braket{pq||rs}+\sum_{p,q}\omega_{pq}(S_{pq}-\delta_{pq})
```
To be able to obtain the gradient, we now must identify the density matrices and then solve the orbital response equations. The density matrices corresponding to the HF contribution are the same as in the previous section. There are additional contributions from the MP2 energy correction, as well as the amplitude response terms. The MP2 energy contribution can be easily identified from the last term of Eq. {eq}`eq:MP2_energy_fct` and gives rise to the following two-particle density matrix: 
```{math}
:label: eq:MP2_2pdm_oovv
\Gamma_{ijab} = -\frac{1}{2} t_{ijab}\,. 
```
The amplitude response ($R^\mathrm{A}_\mathrm{MP2}$) contributions are also reasonably easy to identify:
```{math}
:label: eq:RA_explicit
R^\mathrm{A}_\mathrm{MP2}&=\sum_{i,j,a,b}\tilde{t}_{ijab}\left[t_{ijab}(\epsilon_a+\epsilon_b-\epsilon_i-\epsilon_j)-\braket{ij||ab}\right]\nonumber\\
&=-\sum_{i,j,a,b}\braket{ij||ab}\tilde{t}_{ijab}+\sum_{i,j,a,b}\tilde{t}_{ijab}\left[\sum_{c}\left(\epsilon_{a}\delta_{ac}+\epsilon_b\delta_{bc}\right)t_{ijab}-\sum_k\left(\epsilon_{i}\delta_{ik}+\epsilon_j\delta_{jk}\right)t_{ijab}\right]\nonumber\\
&=-\sum_{i,j,a,b}\braket{ij||ab}\tilde{t}_{ijab}+\sum_{i,j,a,b}\tilde{t}_{ijab}\left(\sum_c f_{ac}t_{ijcb}+\sum_c f_{bc}t_{ijac}-\sum_k f_{ik}t_{kjab}-\sum_k f_{jk}t_{ikab}\right)\nonumber \\
\nonumber \\
&\Downarrow \mathrm{renaming\,\, indices} \nonumber \\
R^\mathrm{A}_\mathrm{MP2}&=-\sum_{i,j,a,b}\braket{ij||ab}\tilde{t}_{ijab}+\sum_{a,b}f_{ab}\sum_{i,j,c}\left(\tilde{t}_{ijac}t_{ijbc}+\tilde{t}_{ijbc}t_{ijac}\right)\nonumber\\
&-\sum_{i,j}f_{ij}\sum_{k,a,b}\left(\tilde{t}_{ikab}t_{jkab}+\tilde{t}_{jkab}t_{ikab}\right)\,,\label{eq:RA_explicit}
```
resulting in the following density matrices:
```{math}
:label: eq:gammaA_ij
\gamma_{ij}^\mathrm{A}=-\sum_{k,a,b}\left(\tilde{t}_{ikab}t_{jkab}+\tilde{t}_{jkab}t_{ikab}\right)\, ,
```
```{math}
:label: eq:gammaA_ab
\gamma_{ab}^\mathrm{A}=\sum_{i,j,c}\left(\tilde{t}_{ijac}t_{ijbc}+\tilde{t}_{ijbc}t_{ijac}\right)\, ,
```
```{math}
:label: eq:GammaA_ijab
\Gamma_{ijab}^\mathrm{A}=-2\,\tilde{t}_{ijab} \, .
```
Combining all density matrices together and replacing the amplitude response Lagrange multipliers with the expression from Eq. {eq}`eq:ampl_rsp_multipliers`, we have:
```{math}
:label: eq:mp2_gamma_ij
\gamma'_{ij}=\gamma_{ij}+\gamma_{ij}^\mathrm{A}=\delta_{ij}-\sum_{k,a,b}\left(\tilde{t}_{ikab}t_{jkab}+\tilde{t}_{jkab}t_{ikab}\right)=\delta_{ij}-\frac{1}{2}\sum_{k,a,b}t_{ikab}t_{jkab} \, , 
```
```{math}
:label: eq:mp2_gamma_ab
\gamma'_{ab}=\gamma_{ab}^\mathrm{A}=\sum_{i,j,c}\left(\tilde{t}_{ijac}t_{ijbc}+\tilde{t}_{ijbc}t_{ijac}\right)=\frac{1}{2}\sum_{i,j,c}t_{ijac}t_{ijbc}\, ,
```
```{math}
:label: eq:mp2_Gamma_ijkl
\Gamma'_{ijkl} = \Gamma_{ijkl} = -2\delta_{ik}\delta_{jl}\, ,
```
```{math}
:label: eq:mp2_Gamma_ijab
\Gamma'_{ijab} = \Gamma_{ijab}+\Gamma_{ijab}^\mathrm{A}=-\frac{1}{2}t_{ijab}-2\,\tilde{t}_{ijab}=-t_{ijab} \, . 
```
Finally, to determine the $\lambda$ orbital response Lagrange multipliers we insert these density matrices into Eq. {eq}`eq:OrbRspEq`. The only non-zero block of $\lambda$ is the occupied-virtual block and the HF density matrices cancel out, so the orbital response equation is:
```{math}
:label: eq:lambda_mp2_ov
&u=i,\, t=a \nonumber\\
&2\lambda_{ia}(\epsilon_i-\epsilon_a)-2\sum_{p,q}(\gamma'_{pq}+\lambda_{pq})\braket{pa||qi}\delta_{i\epsilon_o}\nonumber+\sum_{p,q,r}\left(\Gamma'_{apqr}\braket{ip||qr}-\Gamma'_{ipqr}\braket{ap||qr}\right)=0\nonumber \\
&{^{(1)}}\Leftrightarrow 2\lambda_{ia}(\epsilon_i-\epsilon_a)-2\sum_{j,k}\gamma^\mathrm{A}_{jk}\braket{ja||ki}-2\sum_{b,c}\gamma^\mathrm{A}_{bc}\braket{ba||ci}-2\sum_{j,b}\lambda_{jb}\left(\braket{ja||bi}+\braket{ba||ji}\right)\nonumber\\
&+\sum_{b,j,k}\Gamma'_{abjk}\braket{ib||jk}-\sum_{j,b,c}\Gamma'_{ijbc}\braket{aj||bc}=0\nonumber \\
&{^{(2)}}\Leftrightarrow
\lambda_{ia}(\epsilon_i-\epsilon_a)-\sum_{j,b}\lambda_{jb}\left(\braket{ja||bi}+\braket{ba||ji}\right)=\sum_{j,k}\gamma^\mathrm{A}_{jk}\braket{ki||ja}-\sum_{b,c}\gamma^\mathrm{A}_{bc}\braket{ic||ba}\nonumber \\
&-\frac{1}{2}\sum_{b,j,k}\Gamma'_{abjk}\braket{ib||jk}+\frac{1}{2}\sum_{j,b,c}\Gamma'_{ijbc}\braket{aj||bc}\nonumber \\
&{^{(3)}}\Leftrightarrow \lambda_{ia}(\epsilon_i-\epsilon_a)+\sum_{j,b}\lambda_{jb}\left(\braket{ib||ja}-\braket{ij||ab}\right)= \nonumber \\ 
&=\sum_{j,k}\gamma^{\mathrm{A}}_{jk}\braket{ki||ja}+\sum_{b,c}\gamma^\mathrm{A}_{bc}\braket{ic||ab}-\frac{1}{2}\sum_{b,j,k}\Gamma'_{jkab}\braket{jk||ib}-\frac{1}{2}\sum_{j,b,c}\Gamma'_{ijbc}\braket{ja||bc}
```
Once the $\lambda$ multipliers are determined using an iterative technique, such as the conjugate gradient algorithm, the different blocks of the $\omega$ multipliers can be computed using Eq. {eq}`eq:omega`. Explicitly,
```{math}
:label: eq:omega_mp2_oo
\omega_{ij} = &-\delta_{ij}\epsilon_i - \gamma_{ij}^\mathrm{A}\epsilon_i - \sum_{k,l} \gamma_{kl}^\mathrm{A}\braket{ki||lj} - \sum_{k,a}\lambda_{ka}\left(\braket{ik||ja}+\braket{jk||ia}\right)\nonumber\\
&-\sum_{ab}\gamma_{ab}^\mathrm{A}\braket{ia||jb}-\frac{1}{2}\sum_{k,a,b}\Gamma_{jkab}^\mathrm{A}\braket{ik||ab}\,,\\
```
```{math}
:label: eq:omega_mp2_ov
\omega_{ia} = -\lambda_{ia}\epsilon_i-\frac{1}{2}\sum_{j,k,c}\Gamma^\mathrm{A}_{jkac}\braket{jk||ic} \,,
```
```{math}
:label: eq:omega_mp2_vv
\omega_{ab} = -\gamma_{ab}^\mathrm{A}\epsilon_a-\frac{1}{2}\sum_{i,j,c}\Gamma^\mathrm{A}_{ijbc}\braket{ij||ac}\, .
```

### Excited states

#### Tamm--Dancoff approximation

## Hessians

