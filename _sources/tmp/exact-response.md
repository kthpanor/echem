# Exact states

## Ehrenfest theorem based

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

\begin{equation*}
      \hat{\Omega}_n =  | n \rangle \langle 0 |
\end{equation*}

as well as the corresponding Hermitian conjugate. These operators couple the ground state $|0 \rangle$ to the excited states $| n \rangle$. There will be two equations associated with each excited state, so there are twice as many equations as there are unknown *complex* parameters $P_n(t)$ with independent real and imaginary parts in the time-dependent phase-isolated wave function, $| \bar{\psi}(t) \rangle$.

**3**. Apply perturbation theory:

$$
  P_{n}(t) = P_{n}^{(1)} + P_{n}^{(2)} + P_{n}^{(3)} + \cdots ,
$$

so that the first-order equation will read

$$
\label{eq:pn1}
  \frac{\partial}{\partial t} 
  \langle 0 | [\hat{P}^{(1)} , \hat{\Omega}_n] | 0 \rangle = 
  \frac{1}{i \hbar}
  \langle 0 | [\hat{P}^{(1)} , [\hat{\Omega}_n, \hat{H}_0]] | 0 \rangle 
  - \frac{1}{\hbar}
  \langle 0 | [\hat{\Omega}_n, \hat{V}(t)] | 0 \rangle .
$$

or, equivalently,

$$
  \frac{\partial}{\partial t} P_n^{(1)} =
  -i \omega_{n0} P_n^{(1)} +
  \frac{1}{\hbar} \langle n | \hat{V}(t) | 0 \rangle ,
$$

which, by direct time-integration, yields

\begin{eqnarray*}
\label{first.3}
  P_{n}^{(1)} & = & e^{-i \omega_{n0} t}
  \int^t  \frac{1}{\hbar} \langle n | \hat{V}(t') | 0 \rangle
  e^{i \omega_{n0} t'} dt' \\ \nonumber & = &

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

\begin{eqnarray*}
  \lim_{\omega \rightarrow -\omega_{f0}}
  (\omega_{f0} - \omega) \;
  \langle \langle \hat{\Omega}; \hat{V}^{\omega} \rangle \rangle
  & = &
  \langle 0 | {\hat{\Omega}} | f \rangle
  \langle f | \hat{V}^{\omega} | 0 \rangle .
\end{eqnarray*}

## Complex polarization propagator

### Relaxation in density-matrix theory

The time evolution of the density operator

\begin{equation*}
    \hat{\rho}(t) = | \psi(t) \rangle \langle \psi(t) | =
    | \bar{\psi}(t) \rangle \langle \bar{\psi}(t) | ,
\end{equation*}

is governed by the Liouville equation

\begin{equation*}
\label{liouville}
  \frac{\partial}{\partial t}
  \hat{\rho} = \frac{1}{i\hbar}
  [ \hat{H}, \hat{\rho} ] .
\end{equation*}

The Liouville formalism of quantum mechanics offer the opportunity to include effects of relaxation in the molecular system

\begin{equation*}
\label{damped-liouville}
  \frac{\partial}{\partial t}
  \rho_{mn} = \frac{1}{i\hbar}
  [ \hat{H}, \hat{\rho} ]_{mn} - 
  \gamma_{mn} (\rho_{mn} - \rho_{mn}^\mathrm{eq}) ,
\end{equation*}

where the damping term $\gamma_{mn}$ corresponds to the rate at which matrix element $\rho_{mn}$ of the density operator relaxes to its equilibrium value $\rho_{mn}^\mathrm{eq}$.

### Relaxation in wave-function theory

Expressed in terms of the phase-isolated wave function, the above equation takes the form {cite}`Norman2018`

\begin{align*}
\label{cpp-Ehrenfest}
  \frac{\partial}{\partial t} 
  \langle \bar{\psi} |  \hat{\Omega}_{nm} | \bar{\psi} \rangle = &
  \frac{1}{i\hbar} 
   \langle \bar{\psi} | [\hat{\Omega}_{nm}, \hat{H}] | \bar{\psi} \rangle 
   \\ \nonumber & -
   \gamma_{mn}\left[
     \langle \bar{\psi} |  \hat{\Omega}_{nm} | \bar{\psi} \rangle
   - \langle \bar{\psi}^{\mathrm{eq}} |  \hat{\Omega}_{nm} |
  \bar{\psi}^{\mathrm{eq}} \rangle \right],
\end{align*}

where

\begin{equation*}
      \hat{\Omega}_{nm} =  | n \rangle \langle m | .
\end{equation*}

This equation governs the time-evolution of the wave function with inclusion of relaxation in the system and this formulation of wave mechanics is known as the *complex polarization propagator* (CPP) approach. The CPP equation of motionis solved in the same manner by means of perturbation theory as in the case without damping and the resulting linear response function takes the form

\begin{equation*}
\label{cpp-lrsp}
  \langle \langle \hat{\Omega}; \hat{V}^{\omega} \rangle
  \rangle  = -
  \frac{1}{\hbar} \sum_{n>0} \left[
  \frac{\langle 0 | \hat{\Omega}  | n \rangle 
	\langle n | \hat{V}^{\omega} | 0 \rangle }
       {\omega_{n0}- \omega - i\gamma_{n0}}
+
  \frac{\langle 0 | \hat{V}^{\omega}  | n \rangle 
	\langle n | \hat{\Omega} | 0 \rangle }
       {\omega_{n0} + \omega + i\gamma_{n0}} \right] .
\end{equation*}

(polarizability:label)=
### Polarizability

The polarizability $\alpha_{\alpha\beta}(-\omega; \omega)$ is obtained from the linear response function by letting $\hat{\Omega} = \hat{\mu}_\alpha$ and $\hat{V}^\omega = - \hat{\mu}_\beta$. If we assume a set of common damping parameters ($\gamma_{n0} = \gamma$), the formula for the polarizability takes the form

\begin{eqnarray*}
\label{SOS-alpha}
\alpha_{\alpha\beta} (-\omega;\omega) & = &
\frac{1}{\hbar} \sum_n \left[
\frac{
	\langle 0 | \hat{\mu}_\alpha | n \rangle
	\langle n | \hat{\mu}_\beta | 0 \rangle
      }{
	\omega_{n0} - \omega - i \gamma
        }
+
\frac{
	\langle 0 | \hat{\mu}_\beta | n \rangle
	\langle n | \hat{\mu}_\alpha | 0 \rangle
      }{
	\omega_{n0} + \omega + i \gamma
        } \right] ,
\end{eqnarray*}

or, equivalently, for real wave functions

\begin{equation*}
\label{alpha-complex}
\alpha_{\alpha\beta} (-\omega;\omega) =
\alpha^R_{\alpha\beta} (-\omega;\omega)  + 
i \; \alpha^I_{\alpha\beta} (-\omega;\omega) ,
\end{equation*}

with

\begin{eqnarray*}
\label{alpha-complex-R}
\alpha^R_{\alpha\beta} & = & 
\frac{1}{\hbar} \sum_n \left[
\frac{
	\langle 0 | \hat{\mu}_\alpha | n \rangle
	\langle n | \hat{\mu}_\beta | 0 \rangle (\omega_{n0} - \omega)
      }{
	(\omega_{n0} - \omega)^2 + \gamma^2
        }
\right. \\ \nonumber && + \left.
\frac{
	\langle 0 | \hat{\mu}_\beta | n \rangle
	\langle n | \hat{\mu}_\alpha | 0 \rangle (\omega_{n0} + \omega)
      }{
	(\omega_{n0} + \omega)^2 + \gamma^2
        } \right] ,
\\
\label{alpha-complex-I}
\alpha^I_{\alpha\beta} & = & 
\frac{\gamma}{\hbar} \sum_n \left[
\frac{
	\langle 0 | \hat{\mu}_\alpha | n \rangle
	\langle n | \hat{\mu}_\beta | 0 \rangle
      }{
	(\omega_{n0} - \omega)^2 + \gamma^2
        } 
-
\frac{
	\langle 0 | \hat{\mu}_\beta | n \rangle
	\langle n | \hat{\mu}_\alpha | 0 \rangle
      }{
	(\omega_{n0} + \omega)^2 + \gamma^2
        } \right] .
\end{eqnarray*}

For optical frequencies in the proximity of a transition frequency,
the dispersion of the real and imaginary parts of the polarizability
are dictated by the first terms in the above equations. At least the diagonal components of the $\alpha^R$-tensor are positive below the first resonance. At resonance, with neglect made of the second term, the real part of the polarizability crosses zero and changes sign. The imaginary part of the polarizability, on the other hand, is zero in the limit of static fields, and it has maxima at the resonance frequencies. We also note that in the limit of small frequencies, $\alpha^I$ has a linear dependence on the damping parameter $\gamma$. With respect to sign inversion of the optical frequency, it can be seen that

\begin{eqnarray*}
\alpha^R_{\alpha\beta}(\omega;-\omega) & = &  
\alpha^R_{\alpha\beta}(-\omega;\omega) , \\
\alpha^I_{\alpha\beta}(\omega;-\omega) & = & 
- \alpha^I_{\alpha\beta}(-\omega;\omega) .
\end{eqnarray*}

### Linear absorption cross section

From the imaginary part of the polarizability, we obtain the linear absorption cross section {cite}`Norman2018`

\begin{equation*}\label{cross-section}
\sigma(\omega) = 
\frac{\omega}{\epsilon_0 c}
\mathrm{Im}\left\{
\overline{\alpha}(-\omega;\omega)
\right\} .
\end{equation*}

Complex polarization propagator calculations thus offer a means to determine UV/vis and X-ray absorption spectra by tuning the optical frequency to the region of interest.

