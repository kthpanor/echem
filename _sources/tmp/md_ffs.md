# Force field parameterization

There exist several online tools enabling the generation of force-field parameters for a given molecule, e.g., the Automated Topology Builder (ATB) and Repository [cite ATB]. However, the quality of a generated force field must always be checked with respect to QM calculations. Another way to generate an initial force field is to start from a Gaussian calculation and follow the procedure describe in the Amber tutorial[cite amber] to obtain parameters for the General Amber Force Field (GAFF). You will then need to convert your Amber force field to a Gromacs force field using the ACPYPE Python module [cite actype]. For the present exercise, we will use the standard GAFF force field as well as a refined force field where parameters for bonds, angles, and dihedrals have been modified.

More to be added.
