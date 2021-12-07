# Kohn–Sham method

(kohn-sham-equation)=
## Kohn–Sham equation

Early attempts to employ density functional theory for electronic structure investigations of atoms led to relatively poor results due to inaccurate approximations of the kinetic energy term in the [functional $F$](N-representability). In 1965, Kohn and Sham {cite}`Kohn1965` proposed a solution to this problem by introducing the concept of Kohn–Sham (KS) orbitals and evaluating the kinetic energy according to

\begin{equation*}
T_s[n(\mathbf{r})] = - \frac{\hbar^2}{2 m_\mathrm{e}}
\sum_i \int  \psi_i^*(\mathbf{r}) \nabla^2  \psi_i(\mathbf{r}) \, d^3\mathbf{r}  
\end{equation*}

where the electron density is defined as

\begin{equation*}
n(\mathbf{r}) = \sum_{i} | \psi_i(\mathbf{r})|^2 
\end{equation*}

Underlying this idea, there is a representation of the electronic ground state in terms of a Slater determinant known as the Kohn–Sham reference state

\begin{equation*}
\Psi_s = \frac{1}{\sqrt{N!}} \det \big[\psi_1 \psi_2\cdots\psi_N\big]
\end{equation*}

which ensures that the density is [$N$-representable](N-representability).

To within a small residual correction, the functional $T_s$ amounts to the true kinetic energy functional $T$. Correspondingly, the energy functional can be written

\begin{equation*}
E[n(\mathbf{r})] = T_s[n(\mathbf{r})] + J[n(\mathbf{r})] + E_\mathrm{xc}[n(\mathbf{r})] + \int v(\mathbf{r}) n(\mathbf{r}) \, d^3\mathbf{r}
\end{equation*}

where the electron–electron interaction contribution is partitioned into the classical Coulomb interaction energy

$$
J[n(\mathbf{r})]
= 
\frac{1}{2}
\frac{e^2}{4\pi\varepsilon_0}
\int\int
\frac{n(\mathbf{r}) n(\mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|}
\, d^3\mathbf{r} \, d^3\mathbf{r}'
$$ 

and the exchange–correlation (XC) functional, $E_\mathrm{xc}[n(\mathbf{r})]$, accounting for the non-classical electron–electron interaction and the residual kinetic energy correction, $(T[n(\mathbf{r})] - T_s[n(\mathbf{r})])$.

Within the Kohn–Sham formulation of DFT, the act of minimizing the energy under the constraint of particle conservation closely follows the derivation of the [Hartree–Fock equation](hartree-fock-equation) and the resulting equation from which the best orbitals are found is known as the Kohn–Sham equation. In its canonical form, it reads

\begin{equation*}
\label{eq:ks}
\hat f \psi_i(\mathbf{r}) = \left[ -\frac{\hbar^2}{2 m_\mathrm{e}} \nabla^2 + v_\mathrm{eff}(\mathbf{r}) \right] \psi_i(\mathbf{r}) = \varepsilon_i \psi_i(\mathbf{r})
\end{equation*}

where $\varepsilon_i$ is the KS orbital energy and $v_\mathrm{eff}(\mathbf{r})$ is the effective potential equal to

$$
v_\mathrm{eff}(\mathbf{r}) = 
\frac{e^2}{4\pi\varepsilon_0}
\int \frac{n(\mathbf{r'})}{|\mathbf{r} - \mathbf{r}'|} 
\, d^3\mathbf{r}' 
+ v_\mathrm{xc}(\mathbf{r})
+ v(\mathbf{r})
$$

with

$$
v_\mathrm{xc}(\mathbf{r}) = \frac{\delta E_{xc}}{\delta n(\mathbf{r})}
$$

## Exchange–correlation functionals

The accuracy  of the Kohn-Sham method is determined by the capability of the XC functional to capture the non-classical electron–electron interactions in the molecular system. As the exact form of XC functional is unknown, one introduces physically justified approximations. The development of new XC functionals is an on-going process with new functionals being developed to address the inaccuracies of previous generations.

The XC functionals can be classified by various criteria such as functional parameterization, arguments, etc. Here, we will describe four categories of functionals sorted by argument types and highlight the most prominent representatives of each group. We will assume the XC functional to be written in terms on an *energy density* on the form

$$
E_\mathrm{xc}[n(\mathbf{r})] = \int
n(\mathbf{r}) \,
\varepsilon_\mathrm{xc}(n, \nabla n, \nabla^2 n) \, d^3\mathbf{r} 
$$

where it is implicitly understood that the electron density can be decomposed into $\alpha$- and $\beta$-spin densities that in general may not be equal. 

```{figure} /img/jacobs_ladder.svg
---
name: jacobs_ladder
align: center
---
Jacob's ladder of density functionals.
```

### Local density approximation (LDA)

Functionals of this group depend solely on $n(\mathbf{r})$. In general, they provide a quite poor description of electron-electron interactions in regions of rapidly varying density, like the bonding regions in molecules. Nevertheless, until the late 1980s, almost all DFT calculations of the electronic structure of molecules were carried out using such functionals due to a lack of alternatives. Within LDA, the combination of the Dirac exchange {cite}`fundirac` and the Vosko-Wilk-Nussair correlation {cite}`funvwn` functionals provides the most accurate results in a wide range of applications and it has become recognized as the *de facto* standard LDA functional.

### Generalized gradient approximation (GGA)

Functionals of this group depend on $n(\mathbf{r})$ and $\nabla n(\mathbf{r})$. The inclusion of the electron density gradient enables a more accurate description of electron-electron interactions in regions of rapidly varying density. The first attempts to develop XC functionals of this type were undertaken by Sham {cite}`Sham1971` and Herman and co-workers {cite}`Herman1969` in the 1970s but from a performance point of view they were unsuccessful. Only in 1988, the first reliable exchange functional (B88) of this type was proposed by Becke {cite}`funb88`. Adopting the Becke formulation, several other XC functionals have since been developed. Among these, most notable are the Lee-Yang-Parr (LYP) correlation functional {cite}`funlyp` and the Perdew-Burke-Ernzerhof (PBE) XC functional {cite}`funpbe`. The former LYP functional is an integral component of the most widely used hybrid functional (B3LYP) {cite}`funb3param, funb3lyp`, and the latter PBE functional is the most popular GGA functional in solid-state calculations.  

### Meta generalized gradient approximation (mGGA)

Functionals of this group depend on $n(\mathbf{r})$,  $\nabla n(\mathbf{r})$, and higher-order derivatives of the electron or/and kinetic energy densities. The development of XC functionals of this type has been actively pursued only during the last two decades and after Tao and co-workers introduced the first reliable mGGA functional in 2003 namely Tao-Perdew-Staroverov-Scuseria (TPSS) {cite}`funtpss`.  Among the mGGA XC functionals, the functionals developed by the Truhlar group are the most widely used {cite}`funm06`, and the revised M06-L functional {cite}`funm06l` offers accuracy comparable to hybrid XC functional in many situations.  

### Hybrid exchange-correlation functionals

Functionals of this group explicitly incorporate a fraction of the exact Hartree–Fock exchange in the energy expression. Since the introduction by Becke {cite}`funb3param`, these functionals dominate in calculations of molecular electronic structure. In particular, the  B3LYP functional, which combines Dirac, Becke exchange functionals with VWN, LYP correlation functionals via the three-parameters scheme {cite}`funb3param`, has become the gold standard for DFT calculations of electronic structure and properties of molecules in their ground states. Besides the B3LYP functional, the various hybrid functionals from the Minnesota functional family {cite}`funm06` have reached considerable popularity.

### Range separated exchange-correlation functionals

Functionals of this group explicitly partition the XC contribution to the energy functional into short- and long-range parts. Typically, the same partitioning is also applied to the classical Coulomb term to ensure correct asymptotic behavior of $(J[n(\mathbf{r})] + E_{xc}[n(\mathbf{r})])$. The functionals of this type are frequently employed in studies of excited states and optical properties of molecules in connection with time-dependent DFT methods. Among these functionals, the CAM-B3LYP XC functional {cite}`funcamb3lyp` is the most popular one, and frequently provides a valued alternative to \textit{ab initio} methods in studies of molecular systems with charge-transfer character.  

### Known deficiencies 

After introducing the most popular classes of XC functionals, it is important to mention some of their known deficiencies. 

First, we recognize the inherent inability of standard XC functionals to describe van der Waals interactions. This issue has been addressed in recent years by adding empirical dispersion corrections or by incorporating dispersion corrections via a kinetic-energy density dependent contribution. The former strategy has been adopted by Grimme and co-workers {cite}`grimmed3, grimmed4`, who introduced empirical D3 and D4 models of dispersion interactions and successfully combined these models with several XC functionals. As a result, the B3LYP+D3 XC functional became a new standard for structure optimizations of both weakly and strongly bonded molecular systems. Truhlar and co-workers have adopted the latter strategy to develop the Minnesota functional family {cite}`funm06, funm06l`, and these functionals provide a viable alternative to the D3- or D4-corrected standard functionals.

Second, we draw attention to the fact that there is an amount of self interaction built into the treatment of the electrons as due to the non-exact form of the XC functional. In Hartree–Fock theory, there is a perfect cancellation of the self-interaction terms in between the Coulomb and exchange contributions to the energy, but in DFT this symmetry is broken.