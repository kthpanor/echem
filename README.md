# Compiled homepage
https://kthpanor.github.io/echem/docs/intro.html

# Create a conda environment for the eChem book

Using the echem.yml file (stored at the root of the echem directory) to create a conda environment `echem` will install all needed packages to compile the book.

```
conda env create -f echem.yml
```

# Some commands

This more "manual" installation of specfic packages should no longer be needed.

```
$ pip install -U jupyter-book
$ git clone https://github.com/kthpanor/echem.git
$ cd echem
$ vi docs/dft.md
$ jupyter-book build .
$ open _build/html/index.html
```

**Note**: We need version 0.11 (or higher) to compile our book which is not yet available with `conda install` so therefore do `pip install` as suggested above.

# Publish the html-version

```
$ pip install ghp-import
$ ghp-import -n -p -f _build/html
```

# References
The file `references.bib` in the top directory is a regular BIBTEX file. Add your references in this file. A citation in the text is added with

```
{cite}`Wang2016, Schlegel2011`
```

Multiple lists of references are possible but not yet propoerly implmented.

# Status:
0. Intro
	-  Missing: description of Respondo (AD)
	-  Missing: VIAMD? Gromacs or openmm?
	-  Missing? how to install to run self, yml
	-  Missing? description of structure of echem
	-  To improve: description of multipsi (MD)

1. Tutorials and Workflows: 
    -  UV/vis
    		- TF
		- To improve (not using water?)
    -  Catalysis?
		- Missing
    -  Optical activity
    		- PN
		- To finish
    -  Vibrational
    		- IB, MH
		- To finish
    -  Magnetic
    		-  ZR, OV
		- Missing
    -  X-ray
    		- TF
		- To finish (larger systems)
    -  TF, working on UV/vis and part of MD

2. Electronic Structure Theory:
    - General aspects
    		- PN, MD
		- To finish (spin)
    - Wavefunction
	- General
		- PN, MD
		- To finish
	- Hartree—Fock theory
		- PN
		- To review
	- CI
		- MD
		- To review
	- MP
		- MD?
		- To review
	- MR methods
		- MD
		- To review
    - DFT
    	- PN
	- To finish?

3. Potential Energy Surfaces:
	Section				Responsible	Status	
#-----------------------------------------------------------------------
    - General aspects			IB, MH		To finish
    - Molecular structure optimization	IB, MH		To finish
    - Transition-state theory		AD 		Missing (Ch. 14 overleaf)
    - Interpolation			IB, YM		To finish
    - Conical intersections		AD		Missing (Ch. 15 overleaf)
    - Molecular dynamics		TF, PN		To finish

4. Properties:
	Section				Responsible	Status	
#-----------------------------------------------------------------------
    - Response theory
	- Exact states			PN		To finish
	- CI states			MD		Missing
	- SCF states			PN		To review?
    - ADC				AD, IB		To finish
    - Exciton coupling model		XL		To finish

5. Environment:
	Section				Responsible	Status	
#-----------------------------------------------------------------------
	-  Localized properties		PN		Missing
	-  Polarizable embedding	?		To finish

6. Visualization:
	Section				Responsible	Status	
#-----------------------------------------------------------------------
	- Molecular orbitals		?		To finish
	- Densities			PN		To organize?
	- Transition densities, NTOs	TF		Missing


# Random comments

- Figures in png 300 dpi or svg

