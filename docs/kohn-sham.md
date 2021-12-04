# Kohn–Sham method

## Formulation

Early attempts to employ density functional theory for electronic structure investigations of atoms led to relatively poor results due to inaccurate approximations of the kinetic energy term in the energy functional. In 1965, Kohn and Sham {cite}`Kohn1965` proposed a solution to this problem by introducing the concept of Kohn–Sham (KS) orbitals and evaluating the kinetic energy according to

\begin{equation*}
T_s[n(\mathbf{r})] = - \frac{\hbar^2}{2 m_\mathrm{e}}
\sum_i \int  \psi_i^*(\mathbf{r}) \nabla^2  \psi_i(\mathbf{r}) \, d^3\mathbf{r}  
\end{equation*}

where the electron density is defined as

\begin{equation*}
n(\mathbf{r}) = \sum_{i} | \psi_i(\mathbf{r})|^2 
\end{equation*}

Implying that a single Slater determinant 

\begin{equation*}
\Psi_s = \frac{1}{\sqrt{N}} \det [\psi_1(\mathbf{r}) \psi_2(\mathbf{r})...\psi_N(\mathbf{r})]
\end{equation*}

represents the ground state of the molecular system in the basis of KS orbitals $\{\psi_i(\mathbf{r})\}$. The $T_s[n(\mathbf{r})]$ correspond to the kinetic energy functional within a small residual correction. Taking this into account the energy functional can be rewritten as 

\begin{equation*}
E[n(\mathbf{r})] = T_s[n(\mathbf{r})] + J[n(\mathbf{r})] + E_\mathrm{xc}[n(\mathbf{r})] + \int v(\mathbf{r}) n(\mathbf{r}) d\mathbf{r}
\end{equation*}

where the electron–electron interaction contribution to $E[n(\mathbf{r})]$ is partitioned into the classical Coulomb interaction energy, 
$J[n(\mathbf{r})]$, and exchange–correlation functional, $E_\mathrm{xc}[n(\mathbf{r})]$, which accounts for the non-classical electron–electron interaction and the residual kinetic energy correction. To derive the KS equations, let us consider the minimization of the energy functional 
$E[n(\mathbf{r})]$ with constrain

\begin{equation*}
\int \sum_{i} | \psi_i(\mathbf{r})|^2 d\mathbf{r} = N
\end{equation*}

which implies that the number of electrons should remain constant during the minimization process. To derive the KS method, let us consider the minimization of the energy functional  $E[n(\mathbf{r})]$ with constrain, which implies that the number of electrons should 
remain constant during the minimization process.  Taking the first-order variation of $E[n(\mathbf{r})]$  with respect to KS orbitals and 
equating it zero, one can arrive at the canonical form of the Kohn–Sham equations for the KS orbitals 

\begin{equation*}
\label{eq:ks}
\hat f \psi_i(\mathbf{r}) = \left[ -\frac{1}{2} \nabla^2 + v_\mathrm{eff}(\mathbf{r}) \right] \psi_i(\mathbf{r}) = \varepsilon_i \psi_i(\mathbf{r})
\end{equation*}

where $\varepsilon_i$ is the KS orbital energy and $v_\mathrm{eff}(\mathbf{r})$ is the effective potential composed of the external, Coulomb, and exchange–correlation potentials

\begin{equation*}
    v_\mathrm{eff}(\mathbf{r}) =  v(\mathbf{r}) + \int \frac{\rho(\mathbf{r'})}{| \mathbf{r} - \mathbf{r}'|} d \mathbf{r}' + \frac{\delta E_{xc}[n(\mathbf{r})]}{\delta n(\mathbf{r})} \ . 
\end{equation*}

## Exchange–correlation functionals

The Kohn-Sham method's accuracy is determined by the exchange-correlation functional capability to describe the non-classical electron-electron interaction in a molecular system. Unfortunately, the exact form of exchange-correlation functional is unknown,  and one has to use physically justified approximations. The development of new exchange-correlation functionals is a continuous process, and each year new functionals are developed addressing inaccuracies or deficiencies of the previous generation of exchange-correlation functionals. The exchange-correlation functionals can be classified by various criteria, like functional parameterization, functional arguments, etc. Here, we will describe four groups of exchange-correlation functionals sorted by the type of their arguments and highlight the most prominent ones in each group.

```{figure} /img/jacobs_ladder.svg
---
scale: 80%
name: jacobs_ladder
align: center
---
Jacob's ladder of density functionals.
```

### Local density approximation (LDA)

Functionals of this group are dependent only on the electron density $n(\mathbf{r})$ and provide, in general, a relatively poor description of electron-electron interaction in regions of rapidly varying density, like bonding regions in molecules.  Nevertheless, until the late 1980s, almost all DFT calculations of the electronic structure of molecules have been carried out using such functionals due to a lack of alternatives. Among the various exchange-correlation functionals of this type, the combination of Dirac exchange {cite}{fundirac} and Vosko-Wilk-Nussair correlation {cite}{funvwn} functionals provides the most accurate results in a wide range of applications. It has become recognized as the standard LDA exchange-correlation functional.

### Generalized gradient approximation (GGA)

Functionals of this group are dependent on the electron density and its gradient. The inclusion of electron density gradient dependent contributions in functional form enables a more accurate description of electron-electron interaction in rapidly varying density regions of the molecular system. The first attempt to develop exchange-correlation functionals of this type undertaken by Sham {cite}{Sham1971} and Herman and co-workers {cite}{Herman1969} in the 1970s was unsuccessful. Only in 1988, the first reliable exchange functional B88 of this type was proposed by Becke {cite}{funb88}. Adopting Becke formulation, several new exchange-correlation functionals have been developed in recent decades. Among these functionals, most notable are the Lee-Yang-Parr (LYP) correlation functional {cite}{funlyp} and Perdew-Burke-Ernzerhof (PBE) exchange-correlation functional {cite}{funpbe}. The former LYP functional is an integral component of most popular hybrid exchange-functional B3LYP {cite}{funb3param, funb3lyp}, and the later PBE functional is the most popular GGA functional in solid-state calculations.  

### Meta generalized gradient approximation (mGGA)

Functionals of this group are dependent on $n(\mathbf{r})$,  $\nabla n(\mathbf{r})$, and higher-order derivatives of electron density or/and kinetic energy density. The development of exchange-correlation functionals of this kind has been actively pursued only during the last two decades, after Tao and co-workers introduced the first reliable mGGA functional, namely Tao-Perdew-Staroverov-Scuseria (TPSS), in 2003 {cite}{funtpss}.  Among mGGA exchange-correlation functionals, the functionals developed by the Truhlar group are the most widely used {cite}{funm06}, and the revised M06-L functional {cite}{funm06l} offers accuracy comparable to hybrid exchange-correlation functional in many situations.  

### Hybrid exchange-correlation functionals

Functionals of this group explicitly incorporate a fraction of the exact Hartree-Fock exchange in their energy expression. Since the introduction by Becke {cite}{funb3param}, these functionals dominated the electronic structure calculations of molecules. In particular, the  B3LYP exchange-correlation functional, which combines Dirac, Becke exchange functionals with VWN, LYP correlation functionals via the three-parameters scheme {cite}{funb3param}, has become the de facto golden standard for DFT calculations of electronic structure and properties of molecules in their ground states. Besides B3LYP functional, the various hybrid functionals from the Minnesota functionals family {cite}{funm06} started gaining popularity in recent years.   

### Range separated exchange-correlation functionals

Functionals of this group explicitly partition exchange-correlation contribution to energy functional $E[n(\mathbf{r})]$ into short- and long-range parts. Typically, the same partitioning also applied to Coulomb functional $J[n(\mathbf{r})]$, to ensure correct asymptotic behavior of $J[n(\mathbf{r})] + E_{xc}[n(\mathbf{r})]$. The functionals of this kind are frequently employed in studies of excited states and optical properties of molecules in connection with time-dependent DFT methods. Among these functionals, CAM-B3LYP exchange-correlation functional {cite}{funcamb3lyp} is most popular one, and frequently provides valued alternative to \textit{ab initio} methods in studies of  charge-transfer molecular systems.    

### Known deficiencies 

After introducing the most popular classes of exchange-correlation functionals, we would like to address one deficiency of many exchange-correlation functional, namely the inability to describe van der Waals interactions. This issue has been addressed in recent years by adding empirical dispersion correction to existing exchange-correlation functionals or by incorporating dispersion corrections via kinetic energy density-dependent contribution. The former strategy has been adopted by Grimme and co-workers {cite}{grimmed3, grimmed4}, who introduced empirical D3 and D4 models of dispersion interactions and successfully combined these models with many exchange-correlation functionals. Furthermore, the B3LYP+D3 exchange-correlation functional became a new standard for DFT calculations and is recommended for electronic structure studies of both weakly and strongly bonded molecular systems. Truhlar and co-workers have adopted the latter strategy to develop the Minessota functionals family {cite}{funm06, funm06l}, and these functionals provide a viable alternative to D3 or D4 corrected standard exchange-correlation functionals. 
