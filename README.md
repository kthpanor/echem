# Locaion of the web pages

https://kthpanor.github.io/echem

# Create the echem environment

Using the echem.yml file (stored at the root of the echem directory) to create a conda environment `echem` will install all needed packages to compile the book.

```
conda env create -f echem.yml
```
where the YML file contains

```
name: echem
channels:
  - conda-forge
  - veloxchem
dependencies:
  - python>=3.10
  - jupyter-book
  - jupyterlab
  - jupyterlab-spellchecker
  - jupyterlab_code_formatter
  - webcolors
  - jsonschema-with-format-nongpl
  - black
  - isort
  - ghp-import
  - k3d
  - ipympl
  - ipywidgets
  - openmm
  - py3dmol
  - rdkit
  - veloxchem
  - multipsi
  - dftd4-python
  - xtb-python
  - adcc<0.16
  - pip
  - pip:
    - git+https://github.com/gator-program/respondo.git@v0.0.5
    - git+https://github.com/gator-program/gator.git@vlx-newints
```

# Get a local copy of the book

```
$ git clone 
$ cd echem
$ jupyter-book start
$ open http://localhost:3000/
```
The browser will show the Jupyter book and interactively update it as you edit pages in JupyterLab or any other tool.

# Publish the manual

```
$ git pull
$ git commit -m 'comment on your modifications'
$ git push
```

# Style directives

- Remove input and output of cells with the tags `remove-input` and `remove-output`, respectively.

- Images are included with the syntax:

````
:::{image} ../images/myfig.png
:width: 400px
:align: left
:::
````

- Figures are included with the syntax:

````
:::{figure} ../images/myfig.png
:width: 400px
:align: left

Figure: My figure caption.
:::
````

- Internal links in the manual are created with the syntax:

```
[visible text](#sec:link-name)
```

- Link targets are created with the syntax:

```
(sec:link-name)=
```

# References

The file `references.bib` in the top directory is a regular BIBTEX file. Add your references in this file. A citation in the manual text is added with

```
{cite}`Wang2016, Schlegel2011`
```
