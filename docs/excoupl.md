Exciton coupling model
======================
The exciton model is a practical approach for describing multichromophoric systems.
By dividing the system into individual chromophores, the computational cost is significantly decreased. In the $ab$ $initio$ exciton model {cite}`Li2017`, the Hamiltonian takes the matrix form as in the Frenkel exciton model:

```{math}
:label: eq:exciton_model_matrix
\mathbf{H} = \sum_I^N E_I |\psi_I \rangle \langle \psi_I | +
\sum_{J \neq I}^N V_{IJ} |\psi_I \rangle \langle \psi_J |,
```

where $I$ and $J$ denote the diabatic excited states and $N$ is the total number of states in the system. A diabatic excited state may involve one or more chromophores.

## Subsystems and states

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

where $c_{ia}^A$ is the CI coefficient of the $i \to a$ transition within monomer $A$. In practice $c_{ia}^A$ can be obtained from time-dependent density functional theory calculation with the Tamm$-$Dancoff approximation. Similarly, the wave function for a CT state between monomers $A$ and $B$ can be expressed as:

```{math}
:label: eq:exciton_model_ct_state
|\psi_{CT}^{A \to B} \rangle = |\psi_{h_A l_B} \rangle = \frac{1}{\sqrt{2}} (\{\bar{l_B}^\dagger \bar{h_A}\} + \{l_B^\dagger h_A\}) |\psi_0 \rangle,
```

where $h_A$ denotes the occupied orbital on $A$ and $l_B$ denotes the virtual orbital on $B$.

## Energies and couplings

The excitation energy of the $n$th LE state is expressed as

```{math}
:label: eq:exciton_model_le_energy
E^{A(n)}_{LE} = \sum_{i,a \in A} \sum_{j,b \in A} c_{ia}^A c_{jb}^A & [ \delta_{ij}f_{ab} - \delta_{ab}f_{ij} + 2(ia|jb) - c_{HF}(ij|ab) \\
& + (1-c_{HF})(ia|f_{xc}|jb) ] ,
```

where $c_{HF}$ is the coefficient of Hartree$-$Fock exchange, $f_{xc}$ is the exchange-correlation functional, and $f_{ij}$ and $f_{ab}$ are Fock matrix elements.

Similarly, the excitation energy of a CT state can be expressed as

```{math}
:label: eq:exciton_model_ct_energy
E^{A \to B}_{CT} = & f_{l_Bl_B} - f_{h_Ah_A} + 2(h_Al_B|h_Al_B) - c_{HF}(h_Ah_A|l_Bl_B) \\
& + (1-c_{HF})(h_Al_B|f_{xc}|h_Al_B) .
```

Note that this expression expects that the exchange-correlation functional has correct asymptotic behavior.

We can further write down the expression for the coupling between two LE states on different monomers

```{math}
:label: eq:exciton_model_le_le_coupling
V^{A(m),B(n)}_{LE-LE} = \sum_{i,a \in A} \sum_{j,b \in B} c_{ia}^A c_{jb}^B \left[ 2(ia|jb) - c_{HF}(ij|ab) + (1-c_{HF})(ia|f_{xc}|jb) \right] ,
```

where only two-electron contributions survive.

The coupling between LE state and CT state can be expressed as

```{math}
:label: eq:exciton_model_le_ctAB_coupling
V^{A(n),A \to B}_{LE-CT} = \sum_{i,a \in A} c_{ia}^A & [ \delta_{ih_A}f_{al_B} + 2(ia|h_Al_B) - c_{HF}(ih_A|al_B) \\
& + (1-c_{HF})(ia|f_{xc}|h_Al_B) ] ,
```

and

```{math}
:label: eq:exciton_model_le_ctBA_coupling
V^{A(n),B \to A}_{LE-CT} = \sum_{i,a \in A} c_{ia}^A & [ -\delta_{al_A}f_{ih_B} + 2(ia|h_Bl_A) - c_{HF}(ih_B|al_A) \\
& + (1-c_{HF})(ia|f_{xc}|h_Bl_A) ] ,
```

where the leading contribution comes from one-electron terms (Fock matrix elements).

The coupling between two CT states is

```{math}
:label: eq:exciton_model_ct_ct_coupling
V^{A \to B,B \to A}_{LE-CT} = 2(h_Al_B|h_Bl_A) - c_{HF}(h_Ah_B|l_Bl_A) + (1-c_{HF})(h_Al_B|f_{xc}|h_Bl_A) .
```
