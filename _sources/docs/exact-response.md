Exact states
============
The exact eigenstates of the molecular Hamiltonian are never available in practice and any attempt to derive response functions for exact states may therefore appear as an academic exercise. However, there are several good and practical reasons to do so:

- First and foremost, it reveals the dependence of molecular response properties on intrinsic properties of the system, such as excitation energies and transition moments, and it thereby connects different spectroscopies
- It offers the possibility of constructing few-states-models upon which one can base rational molecular design
- It provides explicit formulas for response functions in the configuration interaction (CI) approximation
- It reveals general properties and symmetries of response functions that are also valid in approximate state theories
- It suggests the identification of excited state properties from a study of poles and residues of *ground state* response functions, e.g., excitation energies are identified as poles of linear response functions
- The formulations of quantum mechanics and perturbation theory that are suitable for approximate state theories are best illustrated in exact state theory

For now we assume the solutions to the eigenvalue problem of $\hat{H}_0$ are known:

$$
  \hat{H}_0 | n \rangle = E_n  | n \rangle ,
$$

where $|n \rangle$ are the exact eigenstates and $E_n$ the respective
energies. Before being exposed to the perturbation, we assume the molecule to be in a reference state $|0\rangle$, in most cases the molecular ground state.

Following the outlined recipe, we find for exact states:

**1**. Parameterization of exact state:
\begin{equation*}
  | \bar{\psi}(t) \rangle = e^{-i\hat{P}(t)} | 0 \rangle ;
\qquad
  \hat{P}(t) = \sum_{n > 0} \left[\rule{0mm}{5mm}
  P_n(t) | n \rangle \langle 0 | +
  P_n^*(t) | 0 \rangle \langle n |
  \right] 
\end{equation*}
%
$$
   e^{-i\hat{P}(t)} |0\rangle = 
   | 0 \rangle \cos \alpha - i \sum_{n > 0} P_n |n\rangle
   \frac{\sin \alpha}{\alpha} ;
   \qquad
   \alpha  = 
   \sqrt{\sum_{n > 0} |P_n|^2} .
$$

**2**. Ehrenfest theorem:
\begin{equation*}
  \label{eq:Ehrenfest}
      \frac{\partial}{\partial t} 
  \langle \bar{\psi}(t) | \hat{\Omega}_n | \bar{\psi}(t) \rangle =
  \frac{1}{i\hbar}
    \langle \bar{\psi}(t) | [\hat{\Omega}_n, \hat{H}] | \bar{\psi}(t)
    \rangle.
\end{equation*}

The time evolution of the exact state will be determined by requiring that this equation is fulfilled for the set of state-transfer operators
%
\begin{equation*}
      \hat{\Omega}_n =  | n \rangle \langle 0 |
\end{equation*}
%
as well as the corresponding Hermitian conjugate. These operators couple the ground state $|0 \rangle$ to the excited states $| n \rangle$. There will be two equations associated with each excited state, so there are twice as many equations as there are unknown *complex* parameters $P_n(t)$ with independent real and imaginary parts in the time-dependent phase-isolated wave function, $| \bar{\psi}(t) \rangle$.

**3**. Apply perturbation theory:

$$
  P_{n}(t) = P_{n}^{(1)} + P_{n}^{(2)} + P_{n}^{(3)} + \cdots ,
$$
%
so that the first-order equation will read
%
$$
\label{eq:pn1}
  \frac{\partial}{\partial t} 
  \langle 0 | [\hat{P}^{(1)} , \hat{\Omega}_n] | 0 \rangle = 
  \frac{1}{i \hbar}
  \langle 0 | [\hat{P}^{(1)} , [\hat{\Omega}_n, \hat{H}_0]] | 0 \rangle 
  - \frac{1}{\hbar}
  \langle 0 | [\hat{\Omega}_n, \hat{V}(t)] | 0 \rangle .
$$
%
or, equivalently,
%
$$
  \frac{\partial}{\partial t} P_n^{(1)} =
  -i \omega_{n0} P_n^{(1)} +
  \frac{1}{\hbar} \langle n | \hat{V}(t) | 0 \rangle ,
$$
%
which, by direct time-integration, yields
%
\begin{eqnarray*}
\label{first.3}
  P_{n}^{(1)} & = & e^{-i \omega_{n0} t}
  \int^t  \frac{1}{\hbar} \langle n | \hat{V}(t') | 0 \rangle
  e^{i \omega_{n0} t'} dt' \\ \nonumber & = &
%
  \frac{1}{i \hbar} 
  \sum_{\omega}
  \frac{\langle n | \hat{V}^{\omega} | 0 \rangle
    F^\omega e^{-i\omega t}}
       {\omega_{n0} - \omega} .
\end{eqnarray*}

**4**. The response functions of an observable $\hat{\Omega}$ are defined by:
\begin{eqnarray*}
\label{eq:rsp}
  \langle \bar{\psi}(t) | \hat{\Omega} | \bar{\psi}(t) \rangle & = &
  \langle 0 | \hat{\Omega} | 0 \rangle + 
  \sum_{\omega}
  \langle \langle \hat{\Omega}; \hat{V}^{\omega} \rangle \rangle 
  F^{\omega}
  e^{-i\omega t}  + \cdots 
\end{eqnarray*}

We identify:

$$
  i \langle 0 | [ \hat{P}^{(1)}, \hat{\Omega}] | 0 \rangle =
  \sum_{\omega}
  \langle \langle \hat{\Omega}; \hat{V}^{\omega} \rangle \rangle 
  F^{\omega}
  e^{-i\omega t} 
$$

or, equivalently,
\begin{equation*}
\label{lrsp}
  \langle \langle \hat{\Omega}; \hat{V}^{\omega} \rangle
  \rangle  = -
  \frac{1}{\hbar} \sum_{n>0} \left[
  \frac{\langle 0 | \hat{\Omega}  | n \rangle 
	\langle n | \hat{V}^{\omega} | 0 \rangle }
       {\omega_{n0}-\omega}
+
  \frac{\langle 0 | \hat{V}^{\omega}  | n \rangle 
	\langle n | \hat{\Omega} | 0 \rangle }
       {\omega_{n0}+\omega} \right] .
\end{equation*}

A residue analysis provides a means to obtain *excited state properties* from a ground state response function. The poles of the linear response function equal excitation energies and the residues are given by
%
\begin{eqnarray*}
  \lim_{\omega \rightarrow -\omega_{f0}}
  (\omega_{f0} - \omega) \;
  \langle \langle \hat{\Omega}; \hat{V}^{\omega} \rangle \rangle
  & = &
  \langle 0 | {\hat{\Omega}} | f \rangle
  \langle f | \hat{V}^{\omega} | 0 \rangle .
\end{eqnarray*}
