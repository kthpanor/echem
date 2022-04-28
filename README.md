<!-- #region -->
# Create a conda environment for the eChem book

Using the echem.yml file (stored at the root of the echem directory) to create a conda environment `echem` will install all needed packages to compile the book.

```
conda env create -f echem.yml
```

Alternatively, you can manually create the environment and install the packages in several steps:

```
conda create -n echem -c conda-forge -c veloxchem python=3.9 veloxchem libblas=*=*mkl
conda activate echem
conda install -c conda-forge jupyter-book matplotlib ghp-import py3dmol h5py k3d
conda install -c conda-forge -c veloxchem multipsi=0.0.1
conda install -c conda-forge -c gator gator
conda install -c conda-forge -c pyscf pyscf
```

# Publish the html-version

```
$ pip install ghp-import
$ ghp-import -n -p -f _build/html
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
        - 9:00-10:00 Introduction, initial questionnaire (experiences, area of expertise, ...) (TF)
        - 10:00-12:00 SCF block (PN, ZR)
        - 13:00-15:00 MP2 block (IB, XL, MH)
        - 15:00-16:00 Wrap-up, block division (TF)
    - 3/5: Theory day
        - 11:00-12:00 Potential energy surface (gradients, Hessians, coordinates, IM) (YM, IB)
        - 13:00-16:00 Group work
        - 16:00-17:00 Wrap-up, presentations
        - 18:00 Dinner at Cyprus(?)
    - 4/5: Tutorial day (Wednesday)
        - 9:00-12:00 Excited states: TDDFT (AD), ADC (MH), exciton coupling (XL, last)
        - 13:00-16:00 Group work
        - 16:00-17:00 Wrap-up, end questionnaire
        
4. Blocks and main responsible
    - Theory
        - MCSCF (MD)
        - Structure optimization (IEB+YMR)
        - Force field parameterization, coordinate transformation, and interpolated PES (IEB+YMR)
        - TDSCF [missing TDDFT (with AD), exercises] (TF+PN)
        - ADC (MH)
    - Tutorials
        - UV/vis (MD+TF)
        - IR/Raman (IEB+MH)
        - X-ray [mainly missing exercises] (MD+TF)

# General status:
0. Intro
	- Missing: description of Respondo (AD)
	- Missing: MD module
    - Missing: Nice main figure with text (IEB)
	- Missing?: Description of structure of echem
	- To improve: description of multipsi (MD)

1. Tutorials and Workflows: 
    -  UV/vis
    	 - TF+MD
    -  Catalysis
         - Add later
    -  Optical activity
         - PN (?)
    -  Vibrational
    	 - IB, MH
    -  Magnetic
    	 - Missing
    -  X-ray
    	 - TF

2. Electronic Structure Theory:
    - General aspects
    	 - PN, MD
    - Wavefunction
	- General
         - PN, MD
	- Hartree—Fock theory
         - PN
	- CI
         - MD
	- MP
         - IEB, MD?
	- MR methods
         - MD
    - DFT
         - PN

3. Potential Energy Surfaces:
    - General aspects
         - IEB, MH
    - Molecular structure optimization
         - IEB, MH
    - Transition-state theory
         - JG, DG
    - Interpolation
         - IEB, YMR
    - Conical intersections
         - AD
    - Molecular dynamics
         - TF

4. Properties:
    - Response theory
	- Exact states
         - PN
	- CI
         - MD
	- SCF
         - PN
    - ADC
         - MH, IEB, (AD)
    - Exciton coupling model
         - XL

5. Environment:
	-  Localized properties
        - PN
	-  Polarizable embedding

6. Visualization:
	- Molecular orbitals
	- Densities
	- Transition densities, NTOs
        - TF
<!-- #endregion -->

# Misc.

- Figures in png 300 dpi or svg

