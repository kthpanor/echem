Exciton coupling model
======================
The exciton model is a practical approach for describing multichromophoric systems.
By dividing the system into individual chromophores, the computational cost is significantly decreased. In the $ab$ $initio$ exciton model, the Hamiltonian takes the matrix form as in the Frenkel exciton model:

```{math}
:label: eq:exciton_model_matrix
\mathbf{H} = \sum_I^N E_I |\psi_I \rangle \langle \psi_I | +
\sum_{J \neq I}^N V_{IJ} |\psi_I \rangle \langle \psi_J |,
```

where $I$ and $J$ denote the diabatic excited states and $N$ is the total number of states in the system. A diabatic excited state may involve one or more chromophores.

## Subsystems and coupling

In a multichromophoric system, it is natural to choose the individual chromophores as subsystems (or monomers). In the $ab$ $initio$ exciton model, several approximations are introduced such that the exciton model Hamiltonian matrix can be constructed from calculations of the chromophore monomers and dimers. These approximations are:

* The monomeric wave functions are approximated by those of the isolated monomers.
* Molecular orbitals from different monomers are approximated to be orthogonal when deriving expressions for couplings.
* Electron repulsion integrals that involve orbitals from three or four monomers are neglected.
* The coupling between the ground state and the excited state is discarded.

The ground state wave function of the system can therefore be constructed as an antisymmetrized product (Slater determinant) of the orbitals of the isolated monomers. Based on the ground state wave function and the above approximations, we can construct two types of diabatic excited states, the locally excited (LE) states and the charge-transfer (CT) excited states. The wave function for a LE state on monomer $A$ can be expressed as:

```{math}
:label: eq:exciton_model_le_state
|\psi_{LE}^{A} \rangle = \sum_{i,a \in A} c_{ia}^{A} |\psi_{ia} \rangle = \frac{1}{\sqrt{2}} \sum_{i,a \in A} c_{ia}^A (\{\bar{a}^\dagger \bar{i}\} + \{a^\dagger i\}) |\psi_0 \rangle,
```

where $c_{ia}^A$ is the CI coefficient of the $i \to a$ transition within monomer $A$. In practice $c_{ia}^A$ can be obtained from time-dependent density functional theory calculation with the Tamm$--$Dancoff approximation. Similarly, the wave function for a CT state between monomers $A$ and $B$ can be expressed as:

```{math}
:label: eq:exciton_model_ct_state
|\psi_{CT}^{A \to B} \rangle = |\psi_{ib} \rangle = \frac{1}{\sqrt{2}} (\{\bar{b}^\dagger \bar{i}\} + \{b^\dagger i\}) |\psi_0 \rangle,
```

where $i$ denotes the occupied orbital on $A$ and $b$ denotes the virtual orbital on $B$.

## States and energies
