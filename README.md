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

