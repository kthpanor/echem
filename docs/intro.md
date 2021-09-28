# Getting started
## Our vision

> For the things we have to learn before we can do them, we learn by doing them.
>
> -Aristotle

Quantum chemistry is a powerful tool.
It is now possible to model complex chemical processes even on a laptop
getting insights into matter at its fundamental scale.

But quantum chemistry is also very complex. Answering a chemical question requires selecting parameters among a wide variety of choices.
Choosing a model system, an electronic structure method, a basis set, a set of properties, and a wide array of parameters which can affect the results in small but sometimes meaningful way... It can be a very daunting task, even for veterans of the field.

Similarly, for those who wish to get a deeper understanding of a method, going through the pages of equation often riddled with inconsistent notations and formulations is very challenging. And at the end, the link between the equation and the computer implementation found in existing softwares can be vague at best.

We believe that a core issue is that humans are not good at learning in abstract terms. We can get very far with a lecture or a textbook, but we will never build as much intuition about how a clock work as by simply breaking one apart and rebuilding it from scratch.

This is exactly the aim of this page, allowing a hands-on approach to computational chemistry.
Together we will dismantle the black box that a computational chemistry code often seems to be, 
go through all the cogs and gears, and build back together some of the main computational methods of modern computational chemistry.
We will do this by presenting the underlying equations, all expressed with consistent notations,
as well as by suggesting a simple python implementation, to really display in action how the theory is implemented into a practical tool.
Additionally, we will put these methods in context by showing how they can be used to address concrete chemical questions,
discussing the strengths and weaknesses of each method and how to best use them to solve practical problems.


## Who we are

The group of people behind the eChem team is engaged in education and research activities in theoretical chemistry. As an integrated component of this work, we develop software for the simulation of complex molecular systems and this e-book lets you explore the fruits of these efforts. We have an overreaching vision and goal to provide:

> A science- and education-enabling software platform for quantum molecular modeling on contemporary and future high-performance computing (HPC) systems
> Software to meet the challenges of the [EuroHPC project](https://eurohpc-ju.europa.eu/)

Behind the term **science-enabling** there are a multitude of software requirements that we find important in our work, including

- coverage of dense 3D system of sizes up to and beyond 500 atoms in the quantum region
- accurate description of electronically excited states that show a more diffuse character than the ground state
- stable and reliable convergence of iterative equation solvers also with use of diffuse basis sets
- time-efficient prototyping of novel scientific approaches
- transparent exposure of data structures to enable in-depth analyses for standard users
- flexible ways to interact with other components of the simulation (such as molecular dynamics, parameterizing the embedding, and data visualization)
- a fast return of results so as to remain in synchronicity with experimental project partners, enabled by an efficient, memory-conservative, hybrid MPI/OpenMP parallelization and SIMD vectorization in core software modules

The term **education-enabling** adds another set of software requirements to this already long list. In this context, the notion of deeper learning refers to taking each student’s understanding of the subject matter to another (deeper) level. Our experience tells us that the process of implementing methods to solve fundamental equations is supremely efficient as a means to achieve that type of deeper learning, but only a small number of students are granted this opportunity as many core modules of scientific software were written a long time ago and have often been made obscure by code optimization. What if we could instead offer access to the needed building blocks to explore quantum chemistry in very much the same manner that we can use the Python NumPy package to explore linear algebra?

## Our software

All the softwares used in these pages are freely available and easy to install, and can serve both as a  development platform,
with user-friendly python interface and easy access to all important variables needed to write a prototype quantum chemistry code,
and as a modern and very efficient code to perform large-scale calculations on high performance computers.

### VeloxChem
x

x

x

x

x

x

x

x


### Gator
x

x

x

x

x

x

x

x


#### adcc
x

x

x

x

x

x

x

x


#### Respondo
x

x

x

x

x

x

x

x


#### HPC-QC
x

x

x

x

x

x

x

x


### Multi-Psi
x

x

x

x

x

x

x

x


### VIAMD
x

x

x

x

x

x

x

x


Other tools
-----------

### Gromacs
x

x

x

x

x

x

x

x


### geomeTRIC
x

x

x

x

x

x

x

x


### PyFrame

The workflow integration between QM and MM/MD

