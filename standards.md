## Make the book as consistent as reasonable

Need to be a bit flexible, as full consistency is not always preferable


- Check units and notation
- Nonzero (not non-zero)
- Ionization energy (not ionization potential)
- Wave function (not wavefunction)
- Use **bold** or *italics* for introducing terms (?)
- Text in equations: $E_{\mathrm{HF}}$
- Generally include visualization of molecule considered, using cell with `remove-input`
- Avoid commas, colon before (and in) equations

Refer to parts of the book as
- part
    - chapter
        - section
            - subsection


### Figures and tables

- Avoid caption if possible (describe in text around instead)
- Static images and Python cells if too expensive
- Figures in png 300 dpi (*potentially* svg)
- Typical sizes
    - matplotlib: 6x4 or 12x12
    - py3Dmol: 400x300
- Typical py3Dmol style
    - `viewer.setViewStyle({"style": "outline", "width": 0.05})`
    - `viewer.setStyle({"stick":{},"sphere": {"scale":0.25}})`


### Standards

```python
scf_results = scf_drv(...)
molecule = vlx.Molecule.from_xyz_string(...)
basis = vlx.MolecularBasis.read(...)
```




### Package abbreviations

```python
import veloxchem as vlx
import numpy as np
import matplotlib.pyplot as plt
import multipsi as mtp
import py3Dmol as p3d
import adcc
import gator
```

### References

Add to .bib (in chronological order) and include URL to create links in the references

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