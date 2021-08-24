<!-- #region -->
# UV/vis absorption/emission

## Spectrum calculation

### Absorption

### Emission

### Vibronic?



## Spectrum analysis



## General aspects

### Comparison to experiment

### Vibrations

### Broadening


## Large systems


## Recommendations







# Old

- Polarization $\rightarrow$ cross-section
- Connection to experiment - Beer-Lambert, 0-0/vertical/adiabatic transition, ...
- Explicit states versus CPP
- Dipole-allowed transitions
- Broadening of spectrum: resolution, life-time, ... Model life-time with Lorentzian/CPP
- Illustration: explicit state versus CPP - small system (water?) and some rather large system
- Analyzing excited states (descriptors...)
- Maybe: comparison of methods and basis set (Table)
- Emission, phosphorescence

Here the overlap of the S<sub>0</sub> vibrational ground state and the third vibrational state of S<sub>i</sub> is such that this transition will contribute more to the intensity than, for instance, vibrational ground state to vibrational ground state. Calculating the different contributions to the different ground states (the Franck–Condon factors) thus yields a good first improvement of the theoretical spectra, featuring a smoother spectrum. If the PES of S<sub>0</sub> and S<sub>i</sub> are sufficiently close, the inclusion of only ground state to ground state may be sufficient, designated 0–0 in the figure. Experimentally, it may only be possible to resolve the 0–0 transition, or the wavelength corresponding to the maximum of absorption, $\lambda$<sub>max</sub> . The vertical transition energy thus gives only a rough estimate of the actual physical situation, but it is in most cases sufficiently accurate.

## Jablonski diagram

## Franck--Condon approximation

## Herzberg--Teller approximation

## Vibronic resolution
<!-- #endregion -->

In this exercise we will consider the five lowest excited states of bi-thiophene molecule in vacuum and water and determine solvent induced changes of excited states using QM/MM methods. 

1. Generate VeloxChem input file for vacuum, and MM potential files (tip3p.pot for TIP3P force field, and pol.pot for polarizable water 
model) using generator.py from GROMACS file. Download generator.py from https://gitlab.com/TheoChemBio/CourseMaterial/mmhpc. 

2. Compute the 5 first excited states for one snapshot in vacuum.

3. Compute the 5 first excited states for one snapshot using TIP3P potential for MM region. Add to @method settings group keyword \textit{potfile: tip3p.pot}

4. Compute the 5 first excited states for one snapshot using polarizable potential for MM region . TIP3P potential for MM region. Add to @method settings group keyword \textit{potfile: pol.pot}

5. Compare solvent induced shifts for all excited states obtained using non-polarizable and polarizable MM regions. 
