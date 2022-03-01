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




### Spacing

- Can do formatting over notebooks
- This one
```python
```

### Naming objects for SCF

- v1
```python
water_xyz = """
C        0.00000000    0.00000000    0.00000000
O        0.00000000    0.00000000    1.43
"""
mol_vlx = vlx.Molecule.read_str(water_xyz) 
bas_vlx = vlx.MolecularBasis.read(mol_vlx,'6-31G') 
```
- v2
```python
mol_str = """
C        0.00000000    0.00000000    0.00000000
O        0.00000000    0.00000000    1.43
"""
molecule = vlx.Molecule.read_str(mol_str, units='angstrom')
basis = vlx.MolecularBasis.read(molecule, 'cc-pVDZ')
```

### Containers

- Standardize SCF part?

- SCF
    - scf_gs = vlx.ScfRestrictedDriver()
    - scf_drv = vlx.ScfRestrictedDriver()
    - scfdrv = vlx.ScfRestrictedDriver()
- Response
    - rsp_drv = vlx.LinearResponseSolver()
    - lrs = vlx.LinearResponseEigenSolver()
- ADC
    - adc_vis = adcc.adc2(scf_gs, n_states=5)
    - adc_state = adcc.adc2(scf_gs, n_states=5)


### Figure generation, size, and showing

- Default colors from plt?

- Same-ish broadening function and parameters?

- Minimal
```python
fig = plt.figure()
```
- More adaptive
```python
plt.figure(figsize = (6,4))
plt.tight_layout(); plt.show()
```
<!-- #endregion -->
