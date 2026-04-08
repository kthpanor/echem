# Locaion of the web pages

https://kthpanor.github.io/echem/

# Environment YML file

To contribute to the manual, create an environment with the Jupyter book software.

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

# Get a local copy of the manual

```
$ git clone https://github.com/kthpanor/echem.git
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

- Remove output from code cells with the tag `remove-output`
- Remove cell input fwith the tag `remove-input`

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

## More cell tags

| Tag | Description |
|----------|-------------|
| remove-cell	| Remove the cell from the rendered output.| 
| remove-input	| Remove the code cell input/source from the rendered output.| 
| remove-output	| Remove the code cell output from the rendered output.| 
| hide-cell	| Hides the cell from the rendered output.| 
| hide-input	| Hides the code cell input/source from the rendered output.| 
| hide-output	| Hides the code cell output from the rendered output.| 
| remove-stderr	| Remove the code cell output stderr from the rendered output. See also project config| 
| remove-stdout	| Remove the code cell output stdout from the rendered output. See also project config| 
| skip-execution	| Skip this cell, when executing the notebook| 
| raises-exception	| Expect the code cell to raise an Exception (and continue execution)| 

# References

The file `references.bib` in the top directory is a regular BIBTEX file. Add your references in this file. A citation in the manual text is added with

```
{cite}`Wang2016, Schlegel2011`
```
