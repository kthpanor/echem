# Create a conda environment for the eChem book

Using the echem.yml file (stored at the root of the echem directory) to create a conda environment `echem` will install all needed packages to compile the book.

```
conda env create -f echem.yml
```

# Commands

Run notebooks by opening *jupyter-lab* (or *jupyter-notebook*):

```
$ jupyter-lab
```

Jupyter can access folders below where it is spawned, so initiate it in, *e.g.* the main eChem folder. Examples on how to use Jupyter can be found [here](https://jupyter-notebook.readthedocs.io/en/latest/examples/Notebook/examples_index.html).


Build the book:

```
$ jupyter-book build .
```

Publish the html-version:

```
$ ghp-import -n -p -f _build/html
```

# Alternative installation

You can manually create the environment and install the packages in several steps:

```
conda create -n echem -c conda-forge -c veloxchem python=3.9 veloxchem libblas=*=*mkl
conda activate echem
conda install -c conda-forge jupyter-book matplotlib ghp-import py3dmol h5py k3d
conda install -c conda-forge -c veloxchem multipsi=0.0.1
conda install -c conda-forge -c gator gator
conda install -c conda-forge -c pyscf pyscf
```

# References
The file `references.bib` in the top directory is a regular BIBTEX file. Add your references in this file. A citation in the text is added with

```
{cite}`Wang2016, Schlegel2011`
```

Multiple lists of references are possible but not yet properly implemented.


# Distribution of labor:

Sections to be added/expanded marked in **bold**.


1. Intro
    - All authors

2. Tutorials and workflows: 
    -  **Photochemistry: MD**
    -  **Vibrational: IEB**
    -  **UV/vis: TF**
    -  **Optical activity: PN**
    -  X-ray: TF

3. Electronic ground states:
    - Primarily PN+MD 
    - **MD will expand CAS parts**

4. Molecular structure and dynamics:
    - Primarily IEB+MH
    - **MD: TF**

5. Spectra and properties:
    - Primarily PN+MD
    - **TDDFT: AD (+TF)**
    - ADC: IEB+MH
    - Exciton coupling: XL

6. Environment:
    - Primarily PN

7. Visualization:
    - Primarily PN+MD+TF