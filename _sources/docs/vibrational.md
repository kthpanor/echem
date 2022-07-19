(sec:vibrational-spec)=
Vibrational spectroscopies
==========================
General aspects on calculating ground-state [infrared (IR)](sec:ir-tutorial) and [Raman](sec:raman-tutorial) spectra of molecules are discussed in this section. The theoretical foundations on which these types of calculations are based are discussed in details in the section on [property gradients](sec:property-gradients),
[Hessians](sec:mol-hessian), and [vibrational analysis](sec:vib-analysis).

```{figure} /img/ir_raman/vibrational_spectroscopy.svg
---
scale: 100%
name: fig-vibspec
align: center
---
Illustration of (a) IR absorption, (b) Rayleigh (elastic) scattering, (c)-(e) Raman scattering: (c) Stokes (the incoming photon energy is larger than the outgoing photon energy) (d) anti-Stokes scattering (the incoming photon energy is smaller than the outgoing photon energy), and (e) resonance Raman scattering. 
```

If the sections on [UV/vis](sec:uv-vis) and [X-ray spectroscopy](sec:x-ray) are mainly focused on transitions between electronic levels, the focus here is on transitions that involve the vibrational fine structure of the molecular spectrum. {numref}`Fig. {number} <fig-vibspec>` illustrates some of the physical processes which can excite molecular vibrational levels. The simplest is the linear absorption of infrared (IR) light which has a photon energy of the same order of magnitude as the energy separation between vibrational levels (see also {numref}`Fig. {number} <fig-radiation>`). Another process that involves the excitation of vibrational degrees of freedom is the inelastic (Raman) scattering of ultraviolet (UV) or visible (vis) monochromatic radiation. Depending on the the relative position of the initial and final vibrational states, Raman scattering is sub-divided into Stokes or anti-Stokes scattering. In Stokes Raman scattering, the initial vibrational state is lower in energy compared to the final state, while in anti-Stokes the initial vibrational state is higher in energy. If the photon energy of the monochromatic radiation is tuned to match the energy of an electronic transition, the inelastic scattering process is called resonance Raman effect.

```{figure} /img/ir_raman/electromagnetic_radiation.svg
---
scale: 100%
name: fig-radiation
align: center
---
The electromagnetic radiation spectrum of interest for IR and Raman spectroscopy (marked in yellow). The inelastic scattering of UV or visible monochromatic radiation is involved in Raman spectroscopy, while the absorption of IR photons is involved in IR absorption spectroscopy. 
```

