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

# Some commands

Build the book:

```
$ jupyter-book build .
```

Publish the html-version:

```
$ ghp-import -n -p -f _build/html
```

# References
The file `references.bib` in the top directory is a regular BIBTEX file. Add your references in this file. A citation in the text is added with

```
{cite}`Wang2016, Schlegel2011`
```

Multiple lists of references are possible but not yet properly implemented.


# Distribution of labor:

0. Intro
    - All authors

1. Tutorials and Workflows: 
    -  UV/vis: TF+MD
    -  Catalysis: Add later
    -  Optical activity: PN (?)
    -  Vibrational: IB+MH
    -  Magnetic: Missing
    -  X-ray: TF

2. Electronic Structure Theory:
    - Primarily PN+MD 

3. Potential Energy Surfaces:
    - Primarily IEB+MH
    - Conical intersections: AD
    - MD: TF

4. Properties:
    - Primarily PN+MD
    - TDDFT: AD+TF
    - ADC: IEB+MH
    - Exciton coupling: XL

5. Environment:
    - Primarily PN

6. Visualization:
    - Primarily PN+MD+TF