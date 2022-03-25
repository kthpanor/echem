Optical activity and dichroism
==============================

The strength in the electronic circular dichroism spectrum associated with transition  $S_n \leftarrow S_0$ is proportional to the rotatory strength

$$
R_{n0} = 
  \sum_{\alpha = x,y,z}
  \mathrm{Im}
  \langle 0 | \hat{\mu}_\alpha | n \rangle
  \langle n | \hat{m}_\alpha | 0 \rangle ,
$$

where $\hat{\mu}$ and $\hat{m}$ are the electric and magnetic dipole operators, respectively.


Continue with the ethylene molecule placed in the $xy$-plane with $z$ as the inter-carbon axis. Perform your calculations at the level of B3LYP/def2-SVP.

1. Stack a second ethylene in a co-facial configuration at a separation distance of 3.5 Å. Perform a standard time-dependent DFT calculation to determine the oscillator and rotatory strengths of Frenkel states $| \sigma \rangle$ and $| \gamma \rangle$. See Figure~1 in Ref.~\cite{Norman2014} for state definitions and compare your excitation energies with those presented in Table~1 in the same reference.

2. Twist the system by 10 degrees and repeat your calculation. Plot the circular dichroism spectrum. Identify the bisignate excitonic band and comment on the sign of the dichroism with respect to your twist orientation. Repeat your calcuation for a twist angle of $-10$ degrees.

3. For the monomer and twisted dimer, trimer, and tetramer, calculate the circular dichroism spectra (1 excited state) using 1, 2, 4, and 8 nodes, respectively. Plot wall times versus number of ethylenes. How does your result compare to the formal scaling of the method?

4. Repeat your calculations with 4 excited states per ethylene, i.e., 4, 8, 12, and 16 excited states for the monomer, dimer, trimer, and tetramer, respectively. How does an increased number of excited states affect wall times? Comment on your results in view of the Davidson numerical algorithm used for finding the lowest eigenvalues and associated eigenvectors.
