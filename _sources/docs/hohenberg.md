# Hohenberg–Kohn theorems

To illustrate the basic idea of DFT, let us consider an isolated molecular system consisting of $A$ nuclei and $N$ electrons. The anti-symmetric $N$-electron wave function describing the electronic ground state of this system, $\Psi_0(\mathbf{r}_1, \mathbf{r}_2,.., \mathbf{r}_N) $, depends on the spatial coordinates of the electrons, $\{\mathbf{r}_i\}$, and it is a solution to the time-independent Schrödinger equation

\begin{equation*}
\label{eq:sch}
\hat{H} \Psi_0 = E_0 \Psi_0 
\end{equation*}

where $E_0$ is the associated electronic energy and $\hat{H}$ is the electronic Hamiltonian

\begin{equation*}
\hat{H} = - \frac{\hbar^2}{2 m_\mathrm{e}}
\sum_i \nabla_i^2 + 
\sum_{i<j}\frac{e^2}{4\pi\varepsilon_0 |\mathbf{r}_i -\mathbf{r}_j|} + 
\sum_i \hat{v}(\mathbf{r}_i)
\end{equation*}

The Hamiltonian is composed of the kinetic energy, electron–electron repulsion, and external potential operators. For an isolated system, the latter reduces to the electron–nuclear attraction operator

\begin{equation*}
\hat{v}(\mathbf{r}_i) = - \sum_K \frac{Z_K e^2}{4\pi\varepsilon_0 
|\mathbf{r}_i - \mathbf{R}_K|}  
\end{equation*}

where $Z_K$ is the proton number of atom $K$ in the system.

Solving the Schrödinger equation for molecular systems is, however, a non-trivial task and requires the introduction of approximations to the $N$-electron wave function. In contrast, the focus in DFT is set on the one-electron density

\begin{equation*}
n_0(\mathbf{r}) = N \int 
|\Psi_0(\mathbf{r}, \mathbf{r}_2,.., \mathbf{r}_N)|^2 
d^3\mathbf{r}_2 \cdots 
d^3\mathbf{r}_N
\end{equation*}

which only depends on three variables namely the Cartesian components of the position vector. 

The first attempts to describe ground states of atoms in terms of electron density were made in the 1920s by Thomas and Fermi. But the use of the density as the primary variable to describe the electronic ground state was not legitimized until 1964, when Hohenberg and Kohn (HK) introduced their theorems for $v$-representable electron densities {cite}`Hohenberg1964`.

```{note}
An electron density is $v$-representable if it is associated with a ground-state solution to the Schrödinger equation.
```

Considering all possible variations in the external potential (system variations), we understand that there is a infinite number of $v$-representable densities associated with $N$-electron systems. We refer to a general density in this set as $n(\mathbf{r})$.

## HK theorem I

> The external potential $\hat{v}(\mathbf{r})$ of a given system is determined to within a trivial additive constant by the $v$-representable electron density of the system.

In other words, there is a unique one-to-one mapping between a ground-state wave function and its one-electron density, and the knowledge of $n_0(\mathbf{r})$ is sufficient to determine the ground-state energy (and other ground-state properties) of the molecular system.

This allows for the separation of the evaluation of the energy into separate the energy functional into two terms

$$
E_0 =
E[n_0(\mathbf{r})] =
F_\mathrm{HK}[n_0(\mathbf{r})] + V[n_0(\mathbf{r})]
$$

where

$$
V[n(\mathbf{r})] =
\int \hat{v}(\mathbf{r}) n(\mathbf{r}) \, d^3\mathbf{r}
$$

and the Hohenberg–Kohn functional is introduced as the sum of kinetic and electron repulsion energies

$$
F_\mathrm{HK}[n(\mathbf{r})] = T[n(\mathbf{r})] + U[n(\mathbf{r})]
$$ 

We note that the HK functional is universal in the sense that it does not depend on the system under study as such dependencies are isolated to the external potential.

## HK theorem II

> For a $v$-representable trial density $n(\mathbf{r})$ fulfilling
>
> \begin{equation*}
\int n(\mathbf{r}) \, d^3\mathbf{r} = N ;
\qquad
n(\mathbf{r}) \geq 0
\end{equation*}
>
> the ground state energy $E_0$ of a molecular system is fulfills and can be detemined from the relation
>
> \begin{equation*}
 E_0 \leq E[n(\mathbf{r})]
\end{equation*}

We recognize this relation as the variational principle in wave function theory. Any practical application of this relation is, however, severely hampered by the fact that the analytical form of the universal HK potential is unknown. Thus, all DFT methods are forced to employ one or another approximate form of it. 

## $N$-representability

The Hohenberg--Kohn theorems provide a theoretical foundation of DFT. Still, they do not give a recipe for the practical implementation of a DFT computational scheme due to the strict requirement of $v$-representability of electron density.  Fortunately, DFT can be reformulated on the grounds of the Hohenberg--Kohn theorems for so-called $N$-representable electron densities.

```{note}
An $N$-representable density is one that satisfies the following conditions 
\begin{equation*}
\int n(\mathbf{r}) \, d^3\mathbf{r} = N ; \quad  
n(\mathbf{r}) \geq 0; \quad 
\int | \nabla n(\mathbf{r})^{1/2}  |^{2} \, d^3\mathbf{r}  \leq \infty
\end{equation*}
```

Employing the minimum energy principle, Levy proved that the universal potential for $N$-representable electron densities could be defined as a constrained search {cite}`Levy1979`

\begin{equation}
F[\rho(\mathbf{r})] = \min_{\Psi \to \rho(\mathbf{r})}  \langle \Psi | \hat T +\hat V_{ee}| \Psi \rangle 
\end{equation}
where $\hat T$ and $\hat V_{ee}$ are the kinetic energy and the electron-electron interaction operators from the Schr{\"o}dinger equation. The functional $F[\rho(\mathbf{r})]$ searches all wave functions $\Psi$ corresponding to the input electron density $\rho(\mathbf{r})$ and leads to minimum of  $\langle \Psi | \hat T +\hat V_{ee}| \Psi \rangle$. Therefore, in the case of $v$-representable electron density, Levy's constrained search allows one to determine the exact universal potential
\begin{equation}
F_\mathrm{HK}[\rho(\mathbf{r})] = F[\rho(\mathbf{r})] \ .
\end{equation}
The Levy constrained allows one to formulate the variational principle for a $N$-representable electron density 
\begin{equation}
E[\rho(\mathbf{r})] = \min_{\rho(\mathbf{r})} (F[\rho(\mathbf{r})] + \int v(\mathbf{r}) \tilde{\rho}(\mathbf{r}) d\mathbf{r}) \ .
\end{equation}
This equation allows to determine the ground state energy of molecular system from the minimum of energy functional $E[\rho(\mathbf{r})]$ 
and provides recipe for obtaining associated $\rho(\mathbf{r})$ via constrained optimization. Introducing a constrain for electron density 
\begin{equation}
\int \rho(\mathbf{r}) d \mathbf{r} = N
\end{equation}
the minimum condition for energy functional $E[\rho(\mathbf{r})]$ can be written as 
\begin{equation}
\delta \{ E[\rho(\mathbf{r})] - \mu \int \rho(\mathbf{r}) d \mathbf{r}\} = 0 \ .     
\end{equation}
Thus, the stationary condition for energy functional $E[\rho(\mathbf{r})]$ can be rewritten as an Euler equation 
\begin{equation}
\mu = v(\mathbf{r}) + \frac{\delta F[\rho(\mathbf{r})]}{\delta \rho(\mathbf{r})} \ , 
\end{equation}
where the Lagrangian multiplier $\mu$ can be interpreted as a chemical potential of the molecular system i.e. $\mu = dE/dN$. 