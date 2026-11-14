# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

The source of the **eChem book**, an online textbook on quantum chemistry and molecular spectroscopy, published at https://kthpanor.github.io/echem/. It is built with Jupyter Book 2 (MyST). Nearly all content lives in Jupyter notebooks (`.ipynb`) that mix theory text (Markdown + LaTeX) with runnable Python examples, mostly using VeloxChem.

## Commands

- Environment: `conda env create -f echem.yml` (or `echem_windows.yml`), then `conda activate echem`. This installs jupyter-book, veloxchem, multipsi, adcc, gator/respondo, xtb, dftd4, openmm, rdkit, py3Dmol and others.
- Live local preview: `jupyter-book start`, then open http://localhost:3000/. It reloads as pages are edited.
- Static build (what CI runs): `jupyter-book build --html`. Output goes to `_build/html`, which is gitignored.

There are no tests and no linter. Pushing to `master` triggers `.github/workflows/` to build and deploy to GitHub Pages. CI does **not** execute notebooks, so the cell outputs saved in each `.ipynb` are exactly what gets published. After you change a code cell, re-run the notebook (for example with `jupyter nbconvert --to notebook --execute --inplace <nb>`) so its stored outputs stay correct, unless the cell is tagged `skip-execution`.

## Structure

- `myst.yml` defines the table of contents and site config. A new page shows up only once you add it to the `toc` here. The book has these parts: `docs/elec_struct` (ground states), `docs/mol_struct` (structure and dynamics), `docs/spec_prop` (spectroscopy and properties), `docs/env` (environment models), `docs/visualize`, and `docs/tutorials`.
- Chapters follow a parent/child pattern. A parent notebook (for example `wft.ipynb` or `xray.ipynb`) is an overview page, and its children are the sections.
- `img/<topic>/` holds static figures. Notebooks reference them with relative paths such as `../../img/misc/foo.svg`.
- `data/<topic>/` holds input and precomputed result files used by notebooks. Note that `.gitignore` excludes `*.h5`, `*.xyz`, `*.pdb`, `*.log`, `*.top`, `*.gro`, `*.itp` and `*.pot`, so new data files of those types must be force-added.
- `docs/references.bib` is the single BibTeX file for the book.
- `standards.md` is the style guide. Its key points are summarized below.

## Authoring conventions

- Cell tags control rendering: `remove-input`, `remove-output`, `remove-cell`, `hide-cell`/`hide-input`/`hide-output`, `skip-execution`, `raises-exception`.
- Expensive calculations appear as static text or as a `skip-execution` cell, with precomputed results loaded from `data/`. A smaller hidden check (`remove-cell`) can verify that the code works.
- Molecule visualizations normally go in a `remove-input` cell using py3Dmol at 400x300, with `viewer.setViewStyle({"style": "outline", "width": 0.05})` and `viewer.setStyle({"stick":{},"sphere": {"scale":0.25}})`.
- Import aliases: `veloxchem as vlx`, `numpy as np`, `matplotlib.pyplot as plt`, `multipsi as mtp`, `py3Dmol as p3d`, `adcc`, `gator`.
- Variable names: `molecule = vlx.Molecule.read_xyz_string(...)`, `basis = vlx.MolecularBasis.read(...)`, `scf_results = scf_drv.compute(...)`.
- MyST syntax: put a link target `(sec:name)=` above a heading and link to it with `[text](#sec:name)`. Cite with ``{cite}`Key1, Key2` ``. Images and figures use `:::{image}` or `:::{figure}` with `:width:` and `:align:`.
- BibTeX entries must include a `url` (DOI link) and match the format of the existing entries.
- Wording: "nonzero", "ionization energy" (not "potential"), "wave function" (two words). Use upright text in math, as in `$E_{\mathrm{HF}}$`. Avoid commas and colons before equations. Introduce new terms in **bold** or *italics*. Prefer explaining a figure in the surrounding text over giving it a caption. Save figures as 300 dpi PNG or as SVG, and make matplotlib plots 6x4 by default.
