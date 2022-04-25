(isr:label)=
## Intermediate state representation

As the name suggest, the intermediate state representation (ISR) approach consists in constructing the ADC matrix with the help of intermediate states $\ket{\tilde{\psi}_I}$. These are obtained by applying excitation operators to the ground state $\ket{0}$. In second quantization, the excitation operator is written as $\hat{C}_I=\{ \hat{a}_a^\dagger\hat{a}_i;\hat{a}_a^\dagger\hat{a}_b^\dagger\hat{a}_j\hat{a}_i, a<b, i<j;... \}$, where the indices $a,b...$ refer to unoccupied orbitals, while $i,j...$ represent occupied orbitals {cite}`Schirmer1991,Mertins1996,Schirmer2004`. Schematic representations of single and double excitations, which are the only two excitation classes that are needed for ADC orders up to ADC(3), are depicted in {numref}`Fig. {number} <fig-isr>`a and {numref}`Fig. {number} <fig-isr>`b, respectively.
```{figure} /img/adc/isr_adc_matrix.svg
---
scale: 100%
name: fig-isr
align: left
---
Illustration of (a) single excitations, (b) double excitations, (c) the structure of the ADC(2) matrix. The numbers in parenthesis indicate the highest order of perturbation theory used to describe each particular block.
```
The intermediate states $\ket{\tilde{\psi}_I}$ are obtained by first applying $\hat{C}_I$ to the many-body ground state:
```{math}
:label: eq:precursor
\ket{\psi_I^{0}}=\hat{C}_I\ket{0} \, , %%%-\ket{0}\bra{0}\hat{C}_I\ket{0} \, ,
```
and then performing a Gram--Schmidt orthogonalization procedure with respect to lower excitation classes (including the ground state)
to obtain precursor states $\ket{\psi_I^{\#}}$, which can then be orthonormalized symmetrically according to {cite}`Wenzel2016`:
```{math}
:label: eq:ISdefinition
\ket{\tilde{\psi}_I}=\sum_J\ket{\psi_J^{\#}}S_{IJ}^{-1/2}\, ,
```
where $S_{IJ}=\braket{\psi_I^{\#}|\psi_J^{\#}}$ are overlap integrals of the precursor states.

The elements of the ADC matrix $\mathbf{M}$ are obtained as matrix elements of the shifted Hamiltonian in the basis of the intermediate states:
```{math}
:label: eq:Mdef
M_{IJ}=\bra{\tilde{\psi}_I}\hat{H}-E_0\ket{\tilde{\psi}_J}=\sum_{K,L}S_{IK}^{-1/2}\bra{ \psi_K^{\#}}\hat{H}-E_0\ket{\psi_L^{\#}}S_{LJ}^{-1/2} \, ,
```
where $E_0$ is the ground state energy.

Having obtained an expression for the ADC matrix, we return to the series expansion of the polarization propagator. In the same way as the propagator is expanded in series, also the matrix elements can be written in terms of orders of perturbation {cite}`Wenzel2016`:
```{math}
:label: eq:Mseries
M_{IJ}^{(k+l+m)}=\sum_{K,L}\left(S_{IK}^{-1/2}\right)^{(k)}\left(\bra{\psi_K^{\#}}\hat{H}-E_0\ket{\psi_L^{\#}}\right)^{(l)}\left(S_{LJ}^{-1/2}\right)^{(m)}\, ,
```
where $k$ and $m$ are the orders of perturbation theory used for the overlap matrices $S_{IK}$ and $S_{JL}$, $l$ is the order used for the matrix elements of the shifted Hamiltonian, and the sum $k+l+m$ represents the order of the contribution to the ADC matrix $\mathbf{M}$.

Using Eq. {eq}`eq:Mseries` in combination with specific classes of excitation operators and truncating the series at the desired order, various levels of ADC theory are obtained. One aspect to note is that the excitation classes needed to construct a specific ADC level are directly connected to the order of perturbation theory. This can be easily seen from {numref}`Fig. {number} <fig-propagator>`, where the zeroth and first order terms are related to single excitations (only one particle-hole pair is involved), while second order terms involve double excitations (two particle-hole pairs are involved). To illustrate this further, we list the explicit expressions for the ADC matrix elements up to second order {cite}`Wenzel2016, Wormit2009`:
%%% %(\textcolor{red}{not completely convinced if to include -- all of -- these})
%%%
```{math}
:label: eq:adcmat_ph_ph_0
M_{ia,jb}^{(0)}=(\epsilon_a-\epsilon_i)\delta_{ab}\delta_{ij}\,,
```
```{math}
:label: eq:adcmat_ph_ph_1
M_{ia,jb}^{(1)}=-\braket{ja||ib}\,,
```
```{math}
:label: eq:adcmat_ph_ph_2
M_{ia,jb}^{(2)}&=\frac{1}{4}\delta_{ij}\sum_{c,k,l}\left[\frac{\braket{ac||kl}\braket{kl||bc}}{\epsilon_a+\epsilon_c-\epsilon_k-\epsilon_l}+\frac{\braket{ac||kl}\braket{kl||bc}}{\epsilon_b+\epsilon_c-\epsilon_k-\epsilon_l}\right]\nonumber\\ &+\frac{1}{4}\delta_{ab}\sum_{c,d,k}\left[\frac{\braket{cd||ik}\braket{jk||cd}}{\epsilon_c+\epsilon_d-\epsilon_i-\epsilon_k}+\frac{{\braket{cd||Ik}\braket{Jk||cd}}}{\epsilon_c+\epsilon_d-\epsilon_j-\epsilon_k}\right] \nonumber \\
&-\frac{1}{2}\sum_{c,k}\left[\frac{{\braket{ac||ik}\braket{jk||bc}}}{\epsilon_a+\epsilon_c-\epsilon_i-\epsilon_k}+\frac{{\braket{ac||ik}\braket{jk||bc}}}{\epsilon_b+\epsilon_c-\epsilon_j-\epsilon_k}\right] \, ,
```
```{math}
:label: eq:adcmat_ph_2p2h_1
M_{ia,klcd}^{(1)}=\braket{kl||id}\delta_{ac}-\braket{kl||ic}\delta_{ad}-\braket{al||cd}\delta_{ik}+{\braket{ak||cd}\delta_{il}}\label{eq:adcmat_ph_2p2h_1}\,,
```
```{math}
:label: eq:c2p2hph_1
M_{ijab,kc}^{(1)}=\braket{kb||ij}\delta_{ac}-\braket{ka||ij}\delta_{bc}-\braket{ab||cj}\delta_{ik}+{\braket{ab||ci}\delta_{jk}}\label{c2p2hph_1}\,,
```
```{math}
:label: eq:adcmat_2p2h_2p2h_0
M_{ijab,klcd}^{(0)}=(\epsilon_a+\epsilon_b-\epsilon_i-\epsilon_j)\delta_{ac}\delta_{bd}\delta_{ik}\delta_{jl}\,,
```
where $\delta_{pq}$ is the Kronecker delta.

The structure of the ADC(2) matrix is depicted in {numref}`Fig. {number} <fig-isr>`c. In principle, the ADC(2) matrix contains all the possible single and double excitations which can be constructed for the system of interest (using a particular basis set). However, except for very small molecules, to include all these excitations would make the ADC matrix intractably large and impossible to diagonalize. In practice, therefore, only the lowest $n$ excited states are ever constructed, where $n$ is the number of states requested by the user. This means that the space of valence excitations is easily accessible, but makes the space of core excitations impossible to reach, except for molecules with very few electrons. An approach to overcome this problem will be discussed in more detail in the next section.

## Core--valence separation
One way to reach the space of core-excitations (and be able to compute, for example, X-ray absorption spectra), without having to deal with an intractably large ADC matrix is provided by the core--valence separation (CVS) approximation. CVS is obtained by decoupling the core and valence excitation spaces and is motivated by the large energy separation between them {cite}`cederbaum1980`. Essentially, it consists of applying only excitation operators that involve one core electron ({numref}`Fig. {number} <fig-cvs>`a and b). This translates into keeping only those blocks of the ADC matrix which include one core orbital ({numref}`Fig. {number} <fig-cvs>`c) and significantly reduces the size of the matrix to be diagonalized ({numref}`Fig. {number} <fig-cvs>`d). The error introduced by this approximation is very small and system-independent {cite}`Herbst2020`.
```{figure} /img/adc/cvs_adc_matrix.svg
---
scale: 100%
name: fig-cvs
align: left
---
Illustration of the (a) single excitations and (b) double excitations involved in the CVS approximation. (c) Schematic representation of the full ADC(2) matrix and (d) the size reduction achieved by the CVS approximation.
```
The CVS approximation allows the computation of X-ray spectroscopies, which are very important techniques for the characterisation of materials, from atoms and molecules, to surfaces and condensed matter systems. In the remaining of this section we will introduce X-ray absorption spectroscopy (XAS) and review some of its features. For a more in-depth discussion see Refs. {cite}`Stohr1992` and {cite}`xrayrev2018`.

The absorption of soft X-ray electromagnetic radiation by a molecular system results in electronic excitation between core initial electronic states (localized, atomic-like) and final states that are delocalized, sometimes continuum-like. A schematic and simplified illustration of this process, alongside the corresponding XAS spectrum is shown in {numref}`Fig. {number} <fig-xas>`. 

The core level binding energies of 1s electrons are element-specific, with very large energy separations between different elements, e.g. $\sim$290 eV (C), $\sim$400 eV (N), and $\sim$500 eV (O). Additionally, atoms of the same species placed in different chemical environments have binding energies which differ by a few electron-volts. This means that absorption peaks corresponding to transitions from chemically inequivalent atoms will occur at different photon energies. The separation between these peaks is called "chemical shift" and it enables the use of XAS for chemical analysis and materials characterization. 
```{figure} /img/adc/xas.svg
---
scale: 100%
name: fig-xas
align: left
---
(a) Schematic potential and (b) schematic representation of the C K-edge X-ray absorption spectrum for the dichloroethylene molecule (adapted from Ref. {cite}`Stohr1992`, Fig. 4.2.).
```
In the case of core-excitations, electron--electron correlation and the presence of the core-hole influence the excitation energies and transition probabilities. These relaxation effects are present also in the case of valence-excitations, but they are more pronounced in the case of core-excitations. The is due the presence of a localized core-hole which leads to a strong net attraction of the electron density towards the probed atom. In addition, the interaction with the excited electron creates a small repulsive polarisation effect in the valence region. Both these two counteracting effects need to be accounted for in order to accurately model XAS. At the ADC level of theory, this is achieved by including at least doubly excited configurations in the construction of the ADC matrix. 

