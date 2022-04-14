X-ray spectroscopies
====================

XAS, XES, XPS, and RIXS in brief. Compare and contrast to UV/vis

- illustration: XAS, XES, XPS, RIXS processes

- illustration: IP of elements / other illustration on element-specificity (?)

- maybe illustration: Augen/fluorescent yield

XAS: focus on NEXAFS (prototypical spectrum, maybe)

In X-ray absorption spectroscopy the photon energy is tuned such that core electrons are targeted and excited to either bound or continuum states, and X-ray emission spectroscopy measures the subsequent decay from such an excited state. These core excitations/de-excitations exhibit strong relaxation effects, making theoretical considerations of the processes particularly challenging. While the removal of a valence electron leaves the remaining electrons relatively unaffected, removing core electrons has a substantial effect on the other electrons due to the significant
change in the screening of the nucleus. Additionally, the core-excited states are embedded in a manifold of valence-excited states that needs to be considered by some computationally feasible method.

```md
struct = gator.get_molecule("""
O       0.0000000000     0.0000000000     0.1187290000
H      -0.7532010000    -0.0000000000    -0.4749160000
H       0.7532010000     0.0000000000    -0.4749160000
""")
basis  = gator.get_molecular_basis(struct,'6-31G')
scfres = gator.run_scf(struct,basis)
adcres = gator.run_adc(struct,basis,scfres,method='cvs-adc2x',singlets=4,core_orbitals=1)
```

or as an input file:

```md
@jobs
task: adc
@end
...
@adc
method: cvs-adc2x
core_orbitals: 1
@end
...
```

[water with 6-311++G**?]

Example input (bash and Python): CVS-ADC(2)-x and CPP-B3LYP and $\Delta$SCF ADC(2) and B3LYP (maybe overlap for B3LYP?)

- illustration: versus experiment (shifted) - XAS and XES

XAS: In continuum
-----------------

CVS and CPP.

- illustration: water spectrum for 0-1000 eV

XES: delta or overlap model
---------------------------

MOM etc

Relaxation effects
------------------

Inclusion through electron correlation, Maybe Z+1. Maybe $\Delta$SCF.

- Illustration: ADC(n) hierarchy (1, 2, 2x, 3 for water)

Discuss difference XAS/XES

TDDFT
-----

## Self-interaction error

While there is no principal differences in considering either the valence or core densities in DFT, the use of approximate functionals introduce an error that is more important to account for when core electrons are addressed: the self-interaction error (SIE). This error is a result of spurious self-interactions for all electrons, but the larger density of the core orbitals makes it more influential for such a considerations — localized densities have larger self-repulsion than delocalized densities. For a single-electron systems the following equality should hold for the two-electron terms:

\begin{equation}
J[\rho] + E_{xc} [\rho] = 0.0
\end{equation}

This cancellation is achieved for, e.g. Hartree–Fock theory, but for any approximate functionals in DFT it will instead result in an erronous contribution to the energy. As the erronous self-interaction occurs in terms of self-repulsion, this will lead to an increase in the energy of the orbital in question (i.e. the absolute value will decrease). The simplest correction to this error was proposed by Perdew and Zunger, by correcting the total exchange-correlation energy by a term

\begin{equation}
XE PZ = E DFT -
(J[\rho_i] + E_{xc}[\rho_i]),
\end{equation}


summed over all singly occupied orbitals. This correction improves the potential energy surfaces, but gives minor improvements for orbital energies. Several schemes of correcting the SIE also exists for correcting orbital energies, 111–113 but none have been utilized in this thesis — we argue instead that the relative energies are most important and thus, in most cases, treat the SIE as a scalar element-dependent error to be accounted for. The SIE will be more thoroughly discussed in Section 5.2, but we here want to reiterate that the self-interaction depends strongly on the density of the orbital, and is thus most significant for core electrons of heavy elements.

DFT-based methods generally suffer from self-interaction error (SIE), which corresponds to the sum of Coulomb and exchange self-interactions that remain because of the use of approximate exchange functionals. In more detail, the Coulomb self-repulsion of each electron included in the Coulomb operator is exactly canceled by the nonlocal exchange self-attraction in Hartree--Fock theory, but this is no longer the case when the exchange operator is replaced by an approximate exchange functional. 128 For occupied orbitals, the SIE is defined as

Imamura and Nakai have determined the ground-state SIE of the C(1s) orbital of CO and shown that it amounts to 0.76, 1.30, and 0.40 eV at the BLYP, 141,142 B3LYP, 143,144 and BHHLYP 145 levels of theory, respectively. 146 The source origin of this SIE in the SCF approximation is the Coulomb operator with its sum over occupied orbitals, and, without a cancellation from the full exact Hartree--Fock exchange operator, it will lead to the SIE given in eq 36. Virtual orbitals will consequently not see an SIE in accordance with that for occupied orbitals. Instead, Imamura and Nakai defined an SIE for virtual orbitals in the context of electronic transitions and which depends on the occupied core orbital from which an excitation takes place.

- Spectrum as function of exact exchange

- maybe table with RPA/TDA

- maybe something with B3LYP versus CAM or something? (Functional dependence)

## Basis set dependence

Z-level, diffuse, tight

- illustration: DZ, TZ, CT, aCT with ADC and B3LYP

## Analysis

Maybe explain input more

Use of specialized CVS space?

Polarization dependence - illustrate

Atom assignment - illustrate

Visualization and descriptors - illustrate/tabulate

CVS versus CPP? Example

## Application? Others. Also check exercises

Maybe: moving to heavier elements: splitting of 2p lines and quadrupole-allowed transitions

Relativistic effects

Maybe RIXS

Consider ethylene, vinylfluoride, difluoroethylene, and trifluoroethylene:

1. Compute the C K-edge X-ray absorption spectra of the ethylene molecule using CVS-ADC(0), CVS-ADC(1), and CVS-ADC(2) in combination with the def2-SVP basis set. Choose an appropriate core space. What are the differences between the spectra computed at the different CVS-ADC levels? Are the differences between levels of theory larger in the case of core excitations or valence excitations? 

2. Compute the C K-edge XAS of vinylfluoride at the CVS-ADC(2) level of theory and using the def2-SVP basis set. Choose an appropriate core space and an appropriate number of excitation vectors. How does your result compare to the calculated spectrum from Ref. \cite{fransson2013}? 
3. Plot your CVS-ADC(2) spectrum in comparison to the experimental spectrum from Ref. \cite{fransson2013} (provided in the file {vinylfluoride.dat}). Identify the contribution of each of the two chemically inequivalent carbon atoms to the experimental peaks.

4. Compute also the C K-edge XAS spectra of difluoroethylene, and trifluoroethylene at the CVS-ADC(2) level of theory. Compare the spectra of ethylene, vinylfluoride, difluoroethylene, and trifluoroethylene to each other. What are the effects of fluorine substitution on the C K-edge spectrum? How large are the chemical shifts?


