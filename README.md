# Create a conda environment for the eChem book

## Use a faster conda solver

With the new conda-libmamba-solver conda can run much faster. Read more in [this blog post](https://www.anaconda.com/blog/conda-is-fast-now).

```
conda info
conda update -n base conda
conda install -n base conda-libmamba-solver
conda config --set solver libmamba
```

## Create the echem environment

Using the echem.yml file (stored at the root of the echem directory) to create a conda environment `echem` will install all needed packages to compile the book.

```
conda env create -f echem.yml
```

## Known issues

- If ``conda update`` does not work as expected ([link to issue](https://github.com/conda/conda/issues/9469)), you can try ``conda install`` with explicit conda version, such as ``conda install -n base conda=23.3.1``

- If you encounter ``InvalidArchiveError`` ([link to issue](https://github.com/conda/conda/issues/12235)), run ``conda clean --all`` and try again.

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

# References
The file `references.bib` in the top directory is a regular BIBTEX file. Add your references in this file (in alphabetic order and with consistent format to existing references). A citation in the text is added with

```
{cite}`Wang2016, Schlegel2011`
```

Multiple lists of references are possible but not yet properly implemented.
