<!-- #region -->
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

Multiple lists of references are possible but not yet properly implemented.


# Course plans:

1. To do
    - Come up with suitable exercises for your section(s)
        - Place exercises at end of subchapters/individual pages
        - PN wants to try that we have solutions just after exercises, for now
    - MD+TF: continue planning course
    - Standardize our Notebooks - see `standards.md`
    - Integrate visualization tools in respective sections
        - Need to get it to run better... Maybe use static images
    - Send out e-mail to potential participants (TF+MD)
    - Meta-lock notebook on week 22, followed by peer review
    - Idea: echem.org

2. General aspects
    - Number of participants (10-15)
    - Expecting participants to install software before
    - Exercises partially run at HPC resources?
    - Questionnaires using Google Forms (?)
    - Might move Wednesday to Thursday (time to reflect, Mickael's schedule)
    - 	Check with Mikael and Markus
    - Worth 1.5 points

3. Schedule: 2/5 - 4/5
    - 2/5: Introduction and theory
        - 9:00-11:00 Introduction, initial questionnaire (experiences, area of expertise, ...)
        - 11:00-12:00 + 13:00-14:00 Hartree-Fock example
        - 14:00-15:00 MP2 exercise
        - 15:00-15:30 Wrap-up, maybe short questionnaire
    - 3/5: Theory day
        - 11:00-12:00 Open questions, potential presentations (start 10:00?)
        - 13:00-15:00 Go over notebooks
        - 15:00-17:00 Exercises
        - 17:00-17:30 Wrap-up, maybe short questionnaire
    - 4/5: Tutorial day
        - 9:00-10:00 Open questions, potential presentations
        - 10:00-12:00 Go over notebooks
        - 13:00-15:00 Exercises
        - 15:00-15:30 Wrap-up, end questionnaire
        
4. Blocks and main responsible
    - Theory
        - DFT (PN+ZR)
        - MCSCF (MD)
        - Structure optimization (IEB+YMR)
        - Force field parameterization and interpolated PES (IEB+YMR)
        - TDSCF (TF) [discuss with PN]
        - ADC (MH)
    - Tutorials
        - UV/vis (MD+TF)
        - IR/Raman (IEB+MH)
        - CD (PN)
        - X-ray (MD+TF)

# Status:
0. Intro
	-  Missing: description of Respondo (AD)
	-  Missing: MD module
    -  Missing: Nice main figure (IEB)
	-  Missing?: Description of structure of echem
	-  To improve: description of multipsi (MD)
	-  To improve: Installation (TF)

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
    	 - ZR, OV
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
    - General aspects
            - IB, MH
            - To finish
    - Molecular structure optimization
            - IB, MH
            - To finish
    - Transition-state theory
            - JG, DG
            - Missing
    - Interpolation
            - IB, YM
            - To finish
    - Conical intersections
            - AD
            - Missing
    - Molecular dynamics
            - TF, PN
            - To finish

4. Properties:
    - Response theory
	- Exact states
            - PN
            - To finish
	- CI states
            - MD
            - Missing
	- SCF states
            - PN
            - To review?
    - ADC
            - AD, IB
            - To finish
    - Exciton coupling model
            - XL
            - To finish

5. Environment:
	-  Localized properties
        - PN
        - Missing
	-  Polarizable embedding
        - To finish

6. Visualization:
	- Molecular orbitals
        - To finish
	- Densities
        - PN
        - To organize?
	- Transition densities, NTOs
        - TF
        - Missing
<!-- #endregion -->

# Random comments

- Figures in png 300 dpi or svg

