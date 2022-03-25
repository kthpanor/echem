Response theory
===============
In a simplistic and generalized view of spectroscopy, the observable is the number of detected particles (e.g. photons or electrons) per unit time in a narrow energy interval and into a small solid angle. The observable is recorded under certain conditions regarding parameters such as temperature, pressure, concentration, and there exists a model (or connection) from which one can deduce molecular properties from the set of measured data. Under typical circumstances in molecular spectroscopies, these connections can be viewed upon as measures of changes in an observable due to the presence of electromagnetic fields with origins attributed to external sources. Compared to atomic fields, the electric fields of conventional lasers are relatively weak. A laser delivering pulses of 10 ns duration and 1 mJ in energy and with a spot size of 100 $\mu$m produces an intensity of about 0.3 GW/cm$^2$. This intensity corresponds to an electric field amplitude of some $F^\omega = 5 \times 10^{-5}$ a.u., which is several orders of magnitude smaller than the internal electric fields that bind electrons in atomic and molecular systems.

Although perturbation theory seems motivated from a comparison of external electric field strengths and internal atomic fields, why do we not simply time propagate the Schrödinger equation for the system in its entire (molecule and time-dependent fields) in a direct and nonperturbational manner? Such a strategy provides the principle for *real-time* methods in quantum chemistry and they have been proven computationally tractable in practical applications, giving direct access to time-dependent observables. In this workshop, we instead insist on invoking perturbation theory and our motivations to do so include:

- The response functions defined in perturbation theory provide the natural meeting point between experiment and theory with a distinct separation of one-, two-, three-photon, etc., optical processes

- Error control is difficult to achieve as the accuracy depends on the propagation scheme, time length, and time step

- It is a numerically elaborate process to separate out nonlinearities from dominant lower-order components in the polarization

- Calculation of vibrational contributions can hardly be made practical in a nonperturbational approach.

Response theory {cite}`Norman2018` may be thought of as a reformulation of standard time-dependent perturbation theory into a form suitable for approximate state theory. Virtually all spectroscopic properties are encompassed by the theory as possible perturbations include:

- time-independent or time-dependent 
- electric or magnetic
- internal or external
- geometric distortions

The derivation of expressions for response functions can appear very different from one source to another and the vast number of technicalities can at first appear overwhelming. There is, however, a pattern to be found and almost like a standard recipe to follow:

1. Find an efficient parameterization of the wave function with careful attention paid to:

	- parameter redundancy
	- state vector normalization
 	- phase isolation and secular divergences

2. Choose an appropriate equation-of-motion that is:

	- based on Schrödinger equation
  	- equivalent in exact state theory
  	- important in approximate state theory

3. Apply perturbation theory to solve for the time-dependent parameters introduced in the parameterization of the wave function. This is where most of the hard work is required but it typically provides limited physical insights

4. Form a well-defined quantity of interest, e.g. the electric dipole moment, $\mu(t)$, and identify response functions in the order expansions of these quantities

In the following, we assume that the system is described by a Hamiltonian of the form

$$
\hat{H} = \hat{H}_0 + \hat{V}(t) ;
\qquad
\hat{V}(t) = \sum_\omega \hat{V}^\omega F^\omega e^{-i\omega t} 
$$

where the summation over optical frequencies includes the positive as well as negative Fourier components of the real periodic external field and $F^\omega$ are the corresponding field amplitudes. The quantum mechanical field-coupling operator $ \hat{V}^\omega$ can in general be frequency dependent but is often not. In the electric-dipole approximation, the coupling between the molecular system and external electric field is provided by minus the electric dipole moment operator, $-\hat{\mu}$.  Before being exposed to the perturbation, we assume the molecule to be in a reference state $|0\rangle$, in most cases the molecular ground state.
