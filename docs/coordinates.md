## Choice of coordinates
A fundamental concept in quantum chemistry is the multi-dimensional potential energy surface (PES). This mathematical object captures the interplay between the electronic and nuclear degrees of freedom in a particular system (molecule, complex, etc.), where the special points on the PES are important in elucidating molecular conformations and explaining the mechanisms of chemical and photo-chemical reactions. The first and second order derivatives of the energy with respect to nuclear displacements (gradients, Hessians, respectively) are key elements in locating these special points (local minima, transition states) and determining minimum energy reaction pathways. In this section, we will discuss how analytical expressions for the first order energy derivatives with respect to nuclear displacements are derived, while in the next section we will see how these energy derivatives are used to perform structure and transition state optimizations.

### Cartesian
There are multiple ways to define relative atomic positions. The simplest way is to use a Cartesian reference system and define each atomic position in terms of its $(x,y,z)$ coordinates. In this coordinate system it is very easy to determine the total energy and energy gradient, but this choice is not favourable for geometry optimization because the energy depends non-linearly on atomic Cartesian coordinates and these coordinates are strongly coupled to each other \cite{Wang2016, Schlegel2011}. A more favourable choice is to work with internal coordinates, such as bond lengths, valence angles and dihedrals (one example is the Z-matrix). Using internal coordinates poses, however, two problems: (1) the choice of coordinates is not unique (they are redundant) and (2) internal coordinates have to be transformed back into Cartesian coordinates to compute the energy and gradient. Several ways of handling these problems and performing geometry optimizations using internal coordinates are discussed below.


 ```{figure} /img/pes/coordinates.svg
:scale: 100%
```

### Redundant coordinates

### Delocalized internal coordinates

### Hybrid delocalized internal coordinates

### Tranlation--rotaion internal coordinates

