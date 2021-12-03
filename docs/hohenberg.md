# Hohenberg–Kohn theorems

To illustrate the basic idea of DFT, let us consider an isolated molecular system with $N$ electrons. The anti-symmetric $N$-electron wave function describing the electronic ground state of this system, $\Psi_0(\mathbf{r}_1, \mathbf{r}_2,.., \mathbf{r}_N) $, depends on the spatial coordinates of the electrons, $\{\mathbf{r}_i\}$, and it is a solution to the time-independent Schrödinger equation

\begin{equation*}
\label{eq:sch}
\hat{H} \Psi_0 = E_0 \Psi_0 
\end{equation*}

where $E_0$ is the associated electronic energy and $\hat{H}$ is the electronic Hamiltonian that is composed of the kinetic energy, electron–electron repulsion, and external potential operators

\begin{equation*}
\hat{H} = \hat{T} + \hat{U} + \hat{V}
\end{equation*}

with

\begin{equation*}
\hat{T} =  - \frac{\hbar^2}{2 m_\mathrm{e}} \sum_i \nabla_i^2 ; \qquad
%
\hat{U} = \sum_{i<j}\frac{e^2}{4\pi\varepsilon_0 |\mathbf{r}_i -\mathbf{r}_j|} ; \qquad
%
\hat{V} = \sum_i \hat{v}(\mathbf{r}_i)
\end{equation*}


For an isolated system, the latter reduces to the electron–nuclear attraction operator

$$
\hat{v}(\mathbf{r}_i) = - \sum_K \frac{Z_K e^2}{4\pi\varepsilon_0 
|\mathbf{r}_i - \mathbf{R}_K|}
$$

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

Considering all possible variations in the external potential (system variations), we understand that there is an infinite number of $v$-representable densities associated with $N$-electron systems. We refer to a general density in this set as $\tilde{n}_0(\mathbf{r})$.

## HK theorem I

> The external potential $\hat{v}(\mathbf{r})$ of a given system is determined to within a trivial additive constant by the $v$-representable electron density of the system.

In other words, there is a unique one-to-one mapping between a ground-state wave function and its one-electron density, and the knowledge of $n_0(\mathbf{r})$ is sufficient to determine the ground-state energy (and other ground-state properties) of the molecular system.

This allows for the separation of the energy into separate the energy functional into two terms {cite}`Hohenberg1964`

$$
E_0 =
E[n_0(\mathbf{r})] =
F_\mathrm{HK}[n_0(\mathbf{r})] + V[n_0(\mathbf{r})]
$$

where

$$
V[n_0(\mathbf{r})] =
\int \hat{v}(\mathbf{r}) n_0(\mathbf{r}) \, d^3\mathbf{r}
$$

and the Hohenberg–Kohn functional is introduced as the sum of kinetic and electron repulsion energies

$$
F_\mathrm{HK}[n_0(\mathbf{r})] = T[n_0(\mathbf{r})] + U[n_0(\mathbf{r})]
$$ 

The HK functional is universal in the sense that it does not depend on the system under study as such dependencies are isolated to the external potential, and it is left undefined for densities that are not $v$-representable.

## HK theorem II

> From the $v$-representable trial densities $\tilde{n}_0(\mathbf{r})$ fulfilling
>
> \begin{equation*}
\int \tilde{n}_0(\mathbf{r}) \, d^3\mathbf{r} = N ;
\qquad
\tilde{n}_0(\mathbf{r}) \geq 0
\end{equation*}
>
> the ground state energy $E_0$ of a molecular system can be detemined from the relation
>
> \begin{equation*}
 E_0 \leq E[\tilde{n}_0(\mathbf{r})]
\end{equation*}

We recognize this relation as the variational principle in wave function theory. Any practical application of this relation in density functional theory is, however, severely hampered by the fact that it is prohibitively difficult to ensure that density variations remain $v$-representable. 

## $N$-representability

The Hohenberg--Kohn theorems provide a theoretical foundation of DFT. Still, they do not give a recipe for the practical implementation of a computational scheme due to the strict requirement of $v$-representability of densities. Fortunately, the theory can be reformulated on the grounds of the Hohenberg--Kohn theorems but for the wider class of so-called $N$-representable densities.

```{note}
An electron density is $N$-representable if it can be obtained from some anti-symmetric wave function. Such a density satisfies the following conditions 
\begin{equation*}
\int n(\mathbf{r}) \, d^3\mathbf{r} = N ; \quad  
n(\mathbf{r}) \geq 0; \quad 
\int \left| \nabla n(\mathbf{r})^{1/2}  \right|^{2} \, d^3\mathbf{r}  \leq \infty
\end{equation*}
```

For a given $N$-representable density, Levy demonstrated that there exist a universal variational functional that delivers the associated sum of kinetic and repulsion energies. This functional can be determined by means of a constrained search over the set of wave functions that yield this density {cite}`Levy1979`

\begin{equation*}
F[n(\mathbf{r})] = \min_{\Psi \to n(\mathbf{r})} \langle \Psi | \hat{T} +\hat{U}| \Psi \rangle 
\end{equation*}

When applied to $v$-representable densities, it was shown that $F$ becomes equal to the HK functional such that

$$
F[n_0(\mathbf{r})] + V[n_0(\mathbf{r})] = E_0
$$

and, for $N$-reprentable densities, it was shown that

$$
F[n(\mathbf{r})] + V[n(\mathbf{r})] \geq E_0
$$

This implies a variational principle with respect to $N$-representable densities and the minimization condition for the energy functional can be written in terms of the variation of a Lagrangian that preserves the number of electrons

\begin{equation*}
\delta \Big[ E[n(\mathbf{r})] + \mu \Big( N - \int n(\mathbf{r}) \, d^3\mathbf{r} \Big) \Big] = 0 
\end{equation*}

where the Lagrange multiplier $\mu$ is the chemical potential of the molecular system, *i.e.*, $\mu = dE/dN$. This stationary condition can alternatively be written

\begin{equation*}
\mu = v(\mathbf{r}) + \frac{\delta F[n(\mathbf{r})]}{\delta n(\mathbf{r})} 
\end{equation*}

