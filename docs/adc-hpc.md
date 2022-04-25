(adc-hpc:label)=
## HPC-QC implementation

To enable large-scale ADC(2) calculations on HPC systems, some considerations need to be taken into account in the implementation. Conventionally, the molecular orbital (MO) integrals, which are essentially tensors with four dimensions, are computed by integral transformation on a single compute node (or workstation). This is unfortunately not practical on HPC clusters, since a cluster node usually has moderate amount of memory which is not enough to store an object with $N^4$ elements (here $N$ is the number of orbitals). Therefore, the MO integral needs to be stored in a distributed way; in other words, each node need to store a portion of the four-dimensional MO integral. It is then natural to distribute the computation of MO integral as well, such that each node only computes part of the integral that needs to be stored locally. This leads to the so-called Fock matrix-driver integral transformation approach, where MO integrals are expressed in terms of Fock-like matrices {cite}`Hohenstein2015, gator`
```{math}
:label: eq:mo-ints-from-fock
\langle ij || kl \rangle = \sum_{\gamma \delta}^{N} \left [ C_{\gamma k} F_{\gamma \delta}^{ij} C_{\delta l} - C_{\gamma l} F_{\gamma \delta}^{ij} C_{\delta k} \right ]
```
where
```{math}
:label: eq:fock-from-dens
F_{\gamma \delta}^{ij} = \sum_{\alpha \beta}^{N} D_{\alpha \beta}^{ij} (\alpha \gamma|\delta \beta)
```
```{math}
:label: eq:dens-from-mo
D_{\alpha \beta}^{ij} = C_{\alpha i} C_{\alpha j}
```
In practice the $ij$ pairs and the corresponding $kl$ matrices are stored on individual nodes. This also facilitates parallelization of subsequent $\sigma$-build (matrix-vector multiplication) with minimal amount of communications.

In HPC-QC implementation of ADC(2), another object that cannot be stored on a cluster node is the excitation vector that involves two particle-hole pairs (namely double excitation). This is because the vector of double excitation has four indices and its sizes increases as $O(N^4)$. Fortunately, explicit storage of the double excitation vector can be avoided by converting the standard eigenproblem in single and double space to a nonlinear eigenproblem in single space {cite}`Hattig2006`. By expanding the ADC(2) equation
```{math}
:label: eq:adc2-mat-vec
\begin{pmatrix} M_{11} & M_{12}\\ M_{21} & M_{22} \end{pmatrix}
\begin{pmatrix} V_{1} \\ V_{2} \end{pmatrix}
= \omega \begin{pmatrix} V_{1} \\ V_{2} \end{pmatrix}
```
we have
```{math}
:label: eq:adc2-mat-vec-v1
M_{11} V_1 + M_{12} V_2 = \omega V_1
```
```{math}
:label: eq:adc2-mat-vec-v2
M_{21} V_1 + M_{22} V_2 = \omega V_2
```
where $\omega$ is the excitation energy, and $V_1$ and $V_2$ are the single and double excitation vectors, respectively. The above two equations can be rewritten as (after eliminating $V_2$)
```{math}
:label: eq:adc2-mat-vec-folded
\left[ M_{11} + M_{12} (\omega - M_{22})^{-1} M_{21} \right] V_1 = \omega V_1
```
Note that $M_{22}$ is diagonal in ADC(2). Let
```{math}
:label: eq:adc2-mat-vec-meff
M^{\rm eff}_{11} = \left[ M_{11} + M_{12} (\omega - M_{22})^{-1} M_{21} \right] V_1
```
we have
```{math}
:label: eq:adc2-mat-vec-nonlinear
M^{\rm eff}_{11} V_1 = \omega V_1
```
This is a nonlinear eigenproblem since $M^{\rm eff}_{11}$ depends on $\omega$.
The excitation energies and single excitation vectors can then be solved
by a nonlinear eigensolver.
Note that the resultant single excitation vectors need to be scaled to satisfy the 
following condition
```{math}
:label: eq:adc2-v1-v2-sum
V_1^{\dagger} V_1^{ } + V_2^{\dagger} V_2^{ } = 1
```
where $V_2^{\dagger} V_2^{ }$ can be calculated as
```{math}
:label: eq:adc2-v2-squared
V_2^{\dagger} V_2^{ } = V_1^{\dagger} M_{12} (\omega - M_{22})^{-2} M_{21} V_1^{ }
```
