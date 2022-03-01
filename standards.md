<!-- #region -->
## Make the book as consistent as reasonable

Need to be a bit flexible, as full consistency is not always preferable

### To do

- Add a spell-checker
- Add any points you think we need to standardize
- Figure out reasonable linking to other parts of book
- Figure out cheap(er) visualizations



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

### Miscellaneous

- Module loads:
    - In the open or in hidden blocks
        - Maybe hidden blocks for tutorials (as not main point there), and open otherwise
- **Bold** or *italics* for introducing terms
- Link between sections
- Static images and Python cells if too expensive

<!-- #endregion -->
