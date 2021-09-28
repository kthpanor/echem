# General aspects

In this section we will cover general aspects of electronic structure theory using wave function theory approaches.

## Hamiltonian
## Variational principle
## Rayleigh--Schrödinger perturbation theory
We want to solve the time-independent Schrödinger equation: 

$$
H \Psi = E \Psi
$$

In perturbation theory, the Hamiltonian is partitioned into *zeroth-order* and *perturbation* terms:

$$
(H_0 + \lambda V) \Psi = E \Psi
$$

with $\lambda$ the coupling strength of the perturbation. We assume that the complete spectrum of the zeroth-order Hamiltonian is known:

$$
\label{E0}
\tag{1}
H_0 \Psi^{(0)}_{n} = E_{n}^{(0)} \Psi^{(0)}_{n}
$$

The energy and wavefunction are expanded in a formal power series in $\lambda$ and terms to the same order are collected on the left- and right-hand sides:

$$
(H_0 + \lambda V) \left( \sum_{k} \lambda^{k} \Psi_{n}^{(k)} \right) = \left( \sum_{k} \lambda^{k} E_{n}^{(k)} \right)  \left(\sum_{j} \lambda^{j} \Psi_{n}^{(j)}\right).
$$

Thus:

$$
(E_{n}^{(0)} - H_{0})\Psi_{n}^{(m)} = V\Psi_{n}^{(m-1)} - \sum_{l = 0}^{m-1} E_{n}^{(m-l)}\Psi_{n}^{(l)}.
$$


Energy corrections and perturbative expansion coefficients, $a_{kn}^{(m)} \equiv \langle \Psi_{k}^{(0)} | \Psi_{n}^{(m)} \rangle$ , for the wavefunction are obtained by projection onto the set of known zeroth-order eigenfunctions:

$$
E_{n}^{(m)} = \left\langle \Psi_{n}^{(0)} \left| V \right| \Psi_{n}^{(m-1)} \right\rangle,
$$

and:

$$
[E_{n}^{(0)} - E_{n}^{(0)}]a_{kn}^{(m)} = 
\sum_{j} \langle \Psi_{k}^{(0)} | V | \Psi_{j}^{(0)} \rangle a_{jn}^{(m-1)} - 
\sum_{l = 0}^{m-1} E_{n}^{(m-l)} a_{kn}^{(l)},\quad 
a_{kn}^{(0)} = \delta_{kn},\,\,a_{nn}^{(m)} = \delta_{m0}
$$

For illustrative purposes, let's show the first terms for state 0 (the ground state). The zeroth order terms $E^{(0)}_{0}$ and $\Psi^{(0)}_{0}$ are defined by equation \ref{E0}. The first order terms are:

$$
E_{0}^{(1)} = \left\langle \Psi_{0}^{(0)} \left| V \right| \Psi_{0}^{(0)} \right\rangle
$$
such that

$$
E_{0}^{(0)} + E_{0}^{(1)} = \left\langle \Psi_{0}^{(0)} \left| H \right| \Psi_{0}^{(0)} \right\rangle
$$

and

$$
\label{PT_amp1}
\tag{2}
a_{k0}^{(1)} = \frac{\langle \Psi_{k}^{(0)} | V | \Psi_{0}^{(0)} \rangle }{E_{k}^{(0)} - E_{0}^{(0)}}
$$

and finally, the second order energy becomes

$$
E_{0}^{(2)} = \left\langle \Psi_{0}^{(0)} \left| V \right| \Psi_{n}^{(1)} \right\rangle
$$

which can be combined with the amplitude equation to become

$$ 
E_{0}^{(2)} =  \sum_k \frac{ \langle \Psi_{0}^{(0)} | V | \Psi_{k}^{(0)} \rangle \langle \Psi_{k}^{(0)} | V | \Psi_{0}^{(0)} \rangle }{E_{k}^{(0)} - E_{0}^{(0)}} 
$$


## Particle densities