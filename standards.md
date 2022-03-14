<!-- #region -->
## Make the book as consistent as reasonable

Need to be a bit flexible, as full consistency is not always preferable

### To do

- Add a spell-checker
- Add any points you think we need to standardize
- Figure out reasonable linking to other parts of book
- Figure out cheap(er) visualizations


### Exercises

- As its own subsection at the end of each section, for example:
    - X-ray spectroscopy
        - Calculating
        - Analysis
        - Exercises

- Solutions on the same page (for starters - might change)


### Equations

- Suppress spin, so $\psi(\mathbf{r},t)$, not $\psi(\mathbf{x},t)$.

- Number equations
    - Prefer not

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

### Object names and code layout

- Largelly follow Patrick's convention (space around `=` and between variables etc)
```
struct = gator.get_molecule(water_xyz)
basis = gator.get_molecular_basis(struct, '6-31G')
scf_drv = vlx.ScfRestrictedDriver()
```
- Has to be quite flexible

### References

Include URL or DOI (but not both) to create links in the references:
```
@article{Sellers1993,
author = {Sellers, Harrell},
title = {The C2-DIIS convergence acceleration algorithm},
journal = {Int. J. Quant. Chem.},
volume = {45},
number = {1},
pages = {31-41},
url = {https://doi.org/10.1002/qua.560450106},
year = {1993}
}
```

or 
```
@article{Sellers1993,
author = {Sellers, Harrell},
title = {The C2-DIIS convergence acceleration algorithm},
journal = {Int. J. Quant. Chem.},
volume = {45},
number = {1},
pages = {31-41},
doi = {10.1002/qua.560450106},
year = {1993}
}
```


### Miscellaneous

- Module loads:
    - In the open or in hidden blocks
        - Maybe hidden blocks for tutorials (as not main point there), and open otherwise
- **Bold** or *italics* for introducing terms
- Link between sections
- Static images and Python cells if too expensive

<!-- #endregion -->
