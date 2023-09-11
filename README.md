# Create a conda environment for the eChem book

## Use a faster conda solver

With the new conda-libmamba-solver conda can run much faster. Read more in [this blog post](https://www.anaconda.com/blog/conda-is-fast-now). For newer distributions libmamba is the default solver, and you can check which one you have and, if needed, update with

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

## Update your echem environment

If the echem.yml file has been updated and you would like to update your echem environment accordingly, use the following command

```
conda env update -f echem.yml --prune
```

## Known issues

- If ``conda update -n base conda`` does not work as expected ([link to issue](https://github.com/conda/conda/issues/9469)), you can try ``conda install`` with explicit conda version, such as ``conda install -n base conda=23.3.1``

- If you encounter ``InvalidArchiveError`` ([link to issue](https://github.com/conda/conda/issues/12235)), run ``conda clean --all`` and try again.

- If you encounter ``DLL load failed`` error on Windows ([link to issue](https://github.com/conda/conda/issues/12161)), try the fix documented in [this link](https://github.com/conda/conda/issues/11795#issuecomment-1335666474).

- When using Python 3.11 you may need to set environment variable ``PYDEVD_DISABLE_FILE_VALIDATION`` to ``1`` (see [this link](https://stackoverflow.com/a/75274358)).

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
