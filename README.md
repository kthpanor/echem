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
0. Intro
	-  IB, description of adcc and geomeTRIC;
	-  AD, description of Respondo;
	-  XL, description of HPC-QC;

1. Tutorials and Workflows: 
    -  TF, working on UV/vis and part of MD;
    -  PN, Optical Activity and Dichroism;
    -  IB, MH finish first draft of Vibrational spectroscopies, fix output and add Raman example;
	-  ZR/OV, add EPR to Magnetic resonances;
	-  X-ray spectroscopies: complete!

2. Electronic Structure Theory:
    -  Hartree—Fock theory: first draft completed;
	-  Configuration Interaction: first draft completed;
	-  M{\o}ller-Plesset: first draft completed;
	-  Multiconfiguration methods: first draft completed;
    -  Density functional theory: PN will add an example of SCF optimization and illustrate the Self-interaction error;

3. Potential Energy Surfaces:
	-  Choice of coordinates: IB will add a notebook example of changing between cartesian and internal coordinates; 
	-  Gradients, Hessians and vibrations: MH will finish first draft;  
    -  Molecular structure optimization: IB will add notebook examples for geometry optimization (XTB, HF, TDHF);
    -  Transition-state theory: AD (Ch. 14 overleaf);
    -  Conical intersections: AD (Ch. 15 overleaf);

4. Molecular dynamics: ML & PN (Ch. 7 overleaf)
	-  YM, IB: potentially add a section on interpolated dynamics, excited state dynamics;

5. Properties:
    - Response theory: MD may add a numerical example for response theory of CI states;
    - Algebraic diagrammatic construction: AD, IB, will add notebook examples for ADC (explicit construction of an ADC matrix, XAS of substituted ethenes?); XL will add a description related to HPC-QC / Large scale calculations;
    - Exciton coupling model: XL will add a notebook example;

6. Environment:
	-  Localized properties: PN will prepare a first draft;
	-  Polarizable embedding:

7. Visualization:
	- Transition densities, NTOs: TF will prepare a first draft.


# Random comments

- Figures in png 300 dpi or svg

- Change gromacs to openmm


# echem 1.0

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
