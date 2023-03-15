## Make the book as consistent as reasonable

Need to be a bit flexible, as full consistency is not always preferable

### To do

- Add any points you think we need to standardize
- Figure out reasonable linking to other parts of book
- Figure out cheap(er) visualizations

### Code formatting

- Install and use the JupyterLab Code Formatter running under `jupyter-lab` but not `jupyter-notebook`

```
conda install -c conda-forge jupyterlab_code_formatter black isort
```

### Spell-checking

- Install and use the JupyterLab Spell-Checker

```
conda install -c conda-forge jupyterlab-spellchecker
```


### Equations

- Use spinors for spin $\psi^\dagger(\mathbf{r})$ and not $\psi^\ast(\mathbf{x})$.

- Avoid numbered equations

- MO indices
    - $ij...$ (occupied)
    - $ab...$ (unoccupied)
    - $pq...$ (general)
    - $tu...$ (active)
- Matrices: $B$ or $\mathbf{B}$
    - Prefer bold 

- Greek letters for atomic orbitals

- Text in equations: $E_{\mathrm{HF}}$

### Terms and abbreviations

- Ionization potential or ionization energy
- Cross section or cross-section
- Wave function or wave-function


### Loading packages (name)

```python
import veloxchem as vlx
import numpy as np
import matplotlib.pyplot as plt
import multipsi as mtp
import py3Dmol as p3d
import adcc
import gator
from pyscf import gto, scf
```

### References

Include URL to create links in the references:
```
@article{Sellers1993,
author  = {Sellers, Harrell},
title   = {The C2-DIIS convergence acceleration algorithm},
journal = {Int. J. Quant. Chem.},
volume  = {45},
pages   = {31-41},
url     = {https://doi.org/10.1002/qua.560450106},
year    = {1993}
}
```


### Miscellaneous

- Module loads:
    - In the open or in hidden blocks
        - Maybe hidden blocks for tutorials (as not main point there), and open otherwise
- **Bold** or *italics* for introducing terms
- Link between sections
- Static images and Python cells if too expensive
- Figures in png 300 dpi (*potentially* svg)
- nonzero, not non-zero
- wave function, not wavefunction

Refer to parts of the book as
- part
    - chapter
        - section
            - subsection