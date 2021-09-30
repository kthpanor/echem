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

# Distribution of labour
0. MD, Getting started intro

1. Tutorials and Workflows: 
    - TF re-structure Xray
    - TF, add first version of UV/vis
    - IB, add vibrational spectroscopies

2. Electronic Structure Theory:
    -  Hartree—Fock theory: (MD proof read)
    -  Density functional theory: PN will move Ch. 5 from overleaf here; will include also a jupyter notebook example for Slater exchange;
    -  Wave-function theory: MD is main responsible and will set up this section (Ch. 6 overleaf) and add “Second quantization”;

3. Potential Energy Surfaces: IB is responsible for this section; will contact ML and MA to ask if they wish/can contribute to the respective subsections;
    -  Molecular structure optimization: IB (lacks examples)
    -  Transition-state theory: AD (Ch. 14 overleaf)
    -  Conical intersections: AD (Ch. 15 overleaf)

4. Molecular dynamics: ML & PN (Ch. 7 overleaf)

5. Properties:
    - Response theory: (OV proof read)
    - Algebraic diagrammatic construction: IB (lacks example)
    - Exciton coupling model: XL (Ch. 11 overleaf)

6. Environment: OV will take main responsibility for this section and will move Ch. 12 from overleaf here, as well as add material on loprop.




#Random comments
 - Figures in png 300 dpi
 - Change gromacs to openmm

#echem 1.0
0. Tutorial
    - UV/vis (TF)
      - pyframe-QM/MM integrated here
    - Catalysis (Mårten? or Oriana?) - need gradients
    - Optical activity and dichroism (PN)
    - Vibrational (MH and IB) - need integral derivatives
    - TM -> OUT
    - magnetic resonances: NMR gone, but EPR yes (OV and ZR)
    - Multi-photon -> OUT
    - X-ray (DONE-ish)
1. Electronic structure
    - Mostly ok
    - Add MCSCF (MD)
    - DFT (PN)
2. PES
    - Conical Intersection -> OUT
    - TST tied to Reactivity
    - MD transfered from overleaf and adapted to openMM (TF and PN)
3. Properties
    - Mostly ok
    - excitation coupling (XL)
4. Environment
    - ESP and RESP : PN will transfer from vlx manual
    - LoProp : PN
    - Polarizable embedding... problematic
5. Visualization
    - Leave as is for moment, maybe integrate to electronic structure?
    - Add somehow orbital (and att/detach densities)vizualisation
