## Choice of coordinates
A fundamental concept in quantum chemistry is the multi-dimensional potential energy surface (PES). This mathematical object captures the interplay between the electronic and nuclear degrees of freedom in a particular system (molecule, complex, etc.), where the special points on the PES are important in elucidating molecular conformations and explaining the mechanisms of chemical and photo-chemical reactions. The first and second order derivatives of the energy with respect to nuclear displacements (gradients, Hessians, respectively) are key elements in locating these special points (local minima, transition states) and determining minimum energy reaction pathways. In this section, we will discuss:
* how the molecular geometry is defined in terms of the atomic positions (link);
* how analytical expressions for the first order energy derivatives with respect to nuclear displacements are derived (link);
* second order derivatives (link);
In the next section, we will see how these energy derivatives are used to perform structure and transition state optimizations (link, link).

### Cartesian
There are multiple ways to define relative atomic positions. The simplest way is to use a Cartesian reference system and define each atomic position in terms of its $(x,y,z)$ coordinates. In this coordinate system it is very easy to determine the total energy and energy gradient, but this choice is not favourable for geometry optimization because the energy depends non-linearly on atomic Cartesian coordinates and these coordinates are strongly coupled to each other {cite}`Wang2016, Schlegel2011`. 

 ```{figure} /img/pes/coordinates.svg
:scale: 100%
```
A more favourable choice is to work with internal coordinates, such as bond lengths, valence angles and dihedrals (one example is the Z-matrix). Using internal coordinates poses, however, two problems: (1) the choice of coordinates is not unique (they are redundant) and (2) internal coordinates have to be transformed back into Cartesian coordinates to compute the energy and gradient. Several ways of handling these problems are discussed below.

### Redundant coordinates
One common choice of coordinates in molecular structure optimization are the Redundant internal coordinates (RIC) {cite}`Pulay1992`.

Given the set of internal coordinates $\mathbf{q}=(q_1,q_2,...,q_{n_q})$, we would like to relate displacements performed in these coordinates to displacements in the $n_x=3N$ Cartesian coordinates $\mathbf{x}=(x_1,x_2,...,x_{n_x})$ (to simplify, we denoted all Cartesian coordinates by $x$). For this purpose, we define a rectangular matrix B (called Wilson B-matrix)  {cite}`orcamanual`:

```{math}
:label: eq:B_ij
    B_{ij} = \frac{\partial q_i}{\partial x_j}\,.
```
%
To determine the changes in $\mathbf{x}$ as a result of changes performed in $\mathbf{q}$, we basically need to obtain the inverse relation $\partial x_j/\partial q_i$. However, because $B$ is rectangular and cannot be directly inverted, we first have to construct a square matrix $\mathbf{G}$:
%
```{math}
:label: eq:G
    \mathbf{G} = \mathbf{B}\mathbf{B}^\mathrm{T} \,.
```
%
Additionally, the internal coordinates are redundant, meaning that the rows of $\mathbf{G}$ are linearly dependent. We therefore have to set up and solve an eigenvalue equation that will allow us to separate out the redundancies {cite}`Pulay1992`:
%
```{math}
	:label: eq:G_eigenval
    \mathbf{G} 
    \begin{pmatrix}
    \mathbf{U} & \mathbf{R}
    \end{pmatrix} = 
    \begin{pmatrix}
    \mathbf{U} & \mathbf{R}
    \end{pmatrix} 
    \begin{bmatrix}
    \boldsymbol{\Lambda} & 0 \\
    0 & 0
    \end{bmatrix} \,.
```
%
This equation has $n=3N-6$ (or $3N-5$ for linear molecules) non-zero eigenvalues and $n_q-n$ zero eigenvalues. The non-zero eigenvalues correspond to linearly independent coordinates, while those which are zero identify the redundant ones. Furthermore, the first $n$ eigenvectors contained in the matrix $\mathbf{U}$ are non-redundant and can be used to define a set of non-redundant coordinates $\mathbf{s}$:
%
```{math}
:label: eq:non-redundant
    \mathbf{s} = \mathbf{U}^\mathrm{T} \mathbf{q} \,. 
```
%
The eigenvalue equation, Eq. {eq}`eq:G_eigenval`, can be used to define a generalized inverse matrix $\mathbf{G}^{-}$:
%
```{math}
:label: eq:generalized_inverse
 \mathbf{G}^{-} = 
    \begin{pmatrix}
    \mathbf{U} & \mathbf{R}
    \end{pmatrix}    
    \begin{bmatrix}
    \boldsymbol{\Lambda}^{-1} & 0 \\
    0\,\,\,\,\,\,\,\, & 0
    \end{bmatrix}
    \begin{pmatrix}
    \mathbf{U}^\mathrm{T}\\ 
    \mathbf{R}^\mathrm{T}
    \end{pmatrix}\, , 
```
%
which, in turn, is used to determine the transformation between Cartesian and internal coordinates displacements {cite}`Wang2016`:
%
```{math}
:label: eq:displacement_q_to_x
    \Delta\mathbf{x} = \mathbf{B}^\mathrm{T}\mathbf{G}^{-}\Delta\mathbf{q}\, . 
```
%
In a similar way, the gradient can be transformed from Cartesian to internal coordinates {cite}`orcamanual`:

```{math}
:label: eq:gradient_x_to_q
    \nabla E_q = \mathbf{G}^{-}\mathbf{B}\nabla E_x \, , 
```

where we have denoted the gradient in Cartesian coordinates by $\nabla E_x$ and the gradient in internal coordinates by $\nabla E_q$. 

With Eqs.{eq}`eq:displacement_q_to_x` and {eq}`eq:gradient_x_to_q` we can now transform a displacement in internal coordinates to a displacement in Cartesian coordinates, compute the energy gradient, and transform the gradient back to internal coordinates.

### Delocalized internal coordinates
Geometry optimizations using internal coordinates require, as we have seen above, to find the generalized inverse matrix $\mathbf{G}$. This is an expensive step, since the matrices involved in Eq. {eq}`eq:generalized_inverse` are constructed in terms of the primitive internal coordinates which are very numerous. One idea to reduce this computational cost is to use the non-redundant coordinates defined in Eq. {eq}`eq:non-redundant`. These new coordinates are actually linear combinations of primitives, and therefore are called "delocalized" internal coordinates (DLC) {cite}`Baker1999`.

Within this coordinate system, we define a new $\mathbf{B}$ matrix which now has a significantly smaller dimension (3$N$-6) compared to the $\mathbf{B}$ matrix in primitive internal coordinates. The new $\mathbf{B}$ matrix is then used to reformulate the displacement and gradient transformations between the Cartesian and delocalized internal coordinate systems. For a detailed derivation of these transformations see Ref. {cite}`Baker1999`.

### Hybrid delocalized internal coordinates
The situation becomes a bit more complicated when the system we would like to optimize is composed of two ore more individual molecules. In this case, fictitious intermolecular bonds, angles, and dihedrals have to be included to define the relative positions between the separate components of the multi-molecular system {cite}`Billeter2000, Wang2016`. In an analogous way as for single molecules, hybrid delocalized internal coordinates (HDLC) can be obtained by including the intermolecular primitives in the transformation defined by Eq. {eq}`eq:non-redundant`. In this approach, the intra- and intermolecular degrees of freedom become mixed, i.e., a displacement in one intermolecular primitive will affects both the intramolecular and intermolecular structures {cite}`Wang2016`.

### Tranlation--rotaion internal coordinates
Another option is to treat the intra- and intermolecular coordinates separately, by introducing translation and rotation coordinates into the set of internal primitives. The resulting choice of coordinates are known as translation--rotation internal coordinates (TRIC) {cite}`Wang2016`. What this actually means is that, instead of adding the full set of Cartesian coordinates in the construction of primitive internal coordinates (HDLC), only six new internal coordinates are used for each molecule to describe its translation and rotation, essentially constraining the intermolecular orientations. This reduces the computational cost and represents a more intuitive way to define the relative orientation between molecules {cite}`Wang2016`.

### References

```{bibliography} ../references.bib
```
