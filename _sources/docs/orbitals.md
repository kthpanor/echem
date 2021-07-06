# Orbitals

## Electronic wave functions

State vectors representing single-electron systems reside in a Hilbert space formed as a direct product

$$
\mathcal{V} = 
\mathcal{V}^\mathrm{o} \otimes 
\mathcal{V}^\mathrm{c} \otimes 
\mathcal{V}^\mathrm{s}
$$

The underlying vector spaces are associated with orbital, charge, and spin degrees of freedom. The elements of $\mathcal{V}$ can be expressed in terms of

$$
\psi^{\mathrm{L}\alpha}(\mathbf{r},t) 
\otimes
\begin{pmatrix}
1\\0
\end{pmatrix}
\otimes
\begin{pmatrix}
1\\0
\end{pmatrix}; \quad
\psi^{\mathrm{L}\beta}(\mathbf{r},t) 
\otimes
\begin{pmatrix}
1\\0
\end{pmatrix}
\otimes
\begin{pmatrix}
0\\1
\end{pmatrix} 
$$

$$
\psi^{\mathrm{S}\alpha}(\mathbf{r},t) 
\otimes
\begin{pmatrix}
0\\1
\end{pmatrix}
\otimes
\begin{pmatrix}
1\\0
\end{pmatrix} ; \quad
\psi^{\mathrm{S}\beta}(\mathbf{r},t) 
\otimes
\begin{pmatrix}
0\\1
\end{pmatrix}
\otimes
\begin{pmatrix}
0\\1
\end{pmatrix} 
$$

with a sum that equals a general wave function for a single-electron system, also known as a spinor,

$$
\psi(\mathbf{r},t) =
\begin{pmatrix}
\psi^{\mathrm{L}\alpha} \\
\psi^{\mathrm{L}\beta} \\
\psi^{\mathrm{S}\alpha} \\
\psi^{\mathrm{S}\beta} 
\end{pmatrix}
$$

In this general case, the particle density becomes

$$
n(\mathbf{r},t)) = 
\psi^\dagger \psi =
\sum_\mathrm{c}^{\mathrm{L,S}}
\sum_\sigma^{\alpha,\beta}
\left|
\psi^{\mathrm{c}\sigma}(\mathbf{r},t) 
\right|^2
$$

- molecular orbitals
	- spin orbitals
	- spatial orbitals
- linear combination of atomic orbitals
