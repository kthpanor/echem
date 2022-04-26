(isr:label)=
## Intermediate state representation


### ISR approach

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
where $E_0$ is the ground-state energy.
This representation of the (shifted) Hamiltonian leads to a Hermitian eigenvalue equation,
```{math}
:label: eq:adc_eigenvalue_eq
  \mathbf{MY} = \mathbf{Y\Omega} \, , \quad \mathbf{Y}^\dagger \mathbf{Y} = \mathbf{1} \, ,
```
the solution of which yields vertical excitation energies ($\omega_n = E_n - E_0$) as eigenvalues collected in the diagonal matrix $\mathbf{\Omega}$,
and the corresponding excitation vectors as eigenvectors $\mathbf{Y}_n$, collected in the columns of $\mathbf{Y}$.

Having obtained an expression for the ADC matrix, we return to the series expansion of the polarization propagator. In the same way as the propagator is expanded in series, also the matrix elements can be written in terms of orders of perturbation {cite}`Wenzel2016`:
```{math}
:label: eq:Mseries
M_{IJ}^{(k+l+m)} = \sum_{K,L}\left(S_{IK}^{-1/2}\right)^{(k)}\left(\bra{\psi_K^{\#}}\hat{H}-E_0\ket{\psi_L^{\#}}\right)^{(l)}\left(S_{LJ}^{-1/2}\right)^{(m)}\, ,
```
where $k$ and $m$ are the orders of perturbation theory used for the overlap matrices $S_{IK}$ and $S_{JL}$,
$l$ is the order used for the matrix elements of the shifted Hamiltonian, and the sum $k + l + m = n$ represents the order of the contribution to the ADC matrix $\mathbf{M}$.
In order to get all contributions of a given order $n$, one needs to sum $k,l,m$ over all terms for which $k+l+m = n$.


### Explicit expressions for ADC(2)

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


### Compactness and separability of the ADC/ISR 


### Comparison to related methods

### Excited-state properties

A distinct advantage of the ISR over the classical propagator approach is that it gives direct access to excited-state wave functions
by expanding it in the intermediate-state basis as
```{math}
:label: eq:isr_wf_expansion
\ket{\Psi_n} = \sum_{J} Y_{Jn} \ket{\tilde{\Psi}_J} \, ,
```
where the elements of the eigenvectors are the expansion coefficients, $Y_{Jn} = \braket{\tilde{\Psi}_J | \Psi_n}$.
This immediately offers the opportunity to calculate physical properties $D_n$ of electronically excited state $n$ via {cite}`Schirmer2004`
```{math}
:label: eq:isr_es_properties
D_n = \bra{\Psi_n} \hat{D} \ket{\Psi_n} = \mathbf{Y}_n^\dagger \, \tilde{\mathbf{D}} \, \mathbf{Y}_n = \sum_{IJ} Y_{In}^* \, \tilde{D}_{IJ} \, Y_{Jn} \, ,
```
where $\tilde{\mathbf{D}}$ is the representation of the operator $\hat{D}$ corresponding to the observable in the intermediate state basis,
```{math}
\tilde{D}_{IJ} = \bra{\tilde{\Psi}_I} \hat{D} \ket{\tilde{\Psi}_J} \, .
```
The matrix $\tilde{\mathbf{D}}$ has a perturbation expansion analogous to that of $\mathbf{M}$ {cite}`Schirmer2004`.
It should be noted that properties calculated in this manner generally differ from those calculated
as derivatives of the energy $E_n$ {cite}`Hodecker2019`.
Transition moments between two different excited states ($m \neq n$) can be obtained in a completely analogous manner as
```{math}
T_{mn} = \bra{\Psi_m} \hat{D} \ket{\Psi_n} = \mathbf{Y}_m^\dagger \, \tilde{\mathbf{D}} \, \mathbf{Y}_n \, .
```
