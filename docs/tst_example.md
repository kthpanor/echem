# Transition state geometry optimization

Once we have a suitable candidate for transition state we can optimize the structure for a transition state. As explained before the transition state will be characterized by a vanishing gradient and one negative Hessian 
eigenvalue. Therefore, in essence, we are exploring the PES looking for a maximum.

Let us try with a simple example:

### Importing modules
``` 
import veloxchem as vlx
```

### Definition of the geometry and setting up the SCF driver
```
geom_str = '''
 C     0.0000000000     0.0000000000     0.0000000000
 N     0.0000000000     0.0000000000     1.1483800000
 H    -1.5853600000     0.0000000000     1.1483800000
'''
geometry = vlx.Molecule.read_str(geom_str, units='angstrom')
basis = vlx.MolecularBasis.read(geometry,"sto-3g")
scf_drv = vlx.ScfRestrictedDriver()
scf_drv.compute(geometry,basis)
```
### Setting up the gradient driver
```
scf_grad_drv = vlx.ScfGradientDriver(scf_drv)
gradient_settings = {'numerical:no'}
scf_grad_drv.update_settings(gradient_settings,method_dict={'xcfun: b3lyp'})
```
### Setting up the optimization driver
Here is where the difference with the geometry optimization to a minimum comes
```
scf_opt_drv = vlx.OptimizationDriver(scf_grad_drv,flag='scf')
scf_opt_drv.transition = True ## Allow the optimization driver to find a TS
scf_opt_drv.hessian = 'first+last' ## Get the hessian for the first and the last geometry so we can get the free energy of the TS
scf_opt_drv.compute(geometry,basis)
```
The output will give four files:
**Geometries during the optimization in xyz format**
```
3    
Iteration 0 Energy -91.47957917
C        0.0000000000    0.0000000000    0.0000000000
N        0.0000000000    0.0000000000    1.1483800000
H       -1.5853600000    0.0000000000    1.1483800000
3    
Iteration 1 Energy -91.48195538
C       -0.0133248403    0.0000000000   -0.0000448664
N        0.0083545166   -0.0000000000    1.1557614176
H       -1.5828748368   -0.0000000000    1.1447117311
3    
Iteration 2 Energy -91.48491343
C       -0.0330446561    0.0000000000    0.0008394403
N        0.0200300505   -0.0000000000    1.1656243432
H       -1.5788526923   -0.0000000000    1.1393381204
3    
Iteration 3 Energy -91.48858856
C       -0.0624148357   -0.0000000000    0.0040155042
N        0.0361560524   -0.0000000000    1.1783350818
H       -1.5722305981   -0.0000000000    1.1313947266

.
.
.

3    
Iteration 15 Energy -91.56485163
C       -0.5126034870   -0.0000000164    0.1647216375
N        0.0565819007   -0.0000000165    1.2461069779
H       -1.3601188371   -0.0000000165    1.0177950052
3    
Iteration 16 Energy -91.56485101
C       -0.5144543278   -0.0000000047    0.1646629279
N        0.0610650638   -0.0000000048    1.2420503655
H       -1.3594737410   -0.0000000048    1.0194003216
3    
Iteration 17 Energy -91.56485103
C       -0.5144063176   -0.0000000047    0.1645795165
N        0.0609018094   -0.0000000048    1.2419920540
H       -1.3594042072   -0.0000000048    1.0194684398
```
**A log file with the output from Geometric**

**The initial and final hessian with the Gibbs free energy**
The last will look like this:
```
# == Summary of harmonic free energy analysis ==
# Note: Rotational symmetry is set to 1 regardless of true symmetry
# Note: Free energy does not include contribution from 1 imaginary mode(s)
# 
# Gibbs free energy contributions calculated at @ 300.00 K:
# Zero-point vibrational energy:                                    7.3990 kcal/mol 
# H   (Trans + Rot + Vib = Tot):   1.4904 +   0.8942 +   0.0003 =   2.3849 kcal/mol 
# S   (Trans + Rot + Vib = Tot):  35.8748 +  16.4958 +   0.0009 =  52.3716 cal/mol/K
# TS  (Trans + Rot + Vib = Tot):  10.7624 +   4.9487 +   0.0003 =  15.7115 kcal/mol 
# 
# Ground State Electronic Energy    : E0                        =   -91.56485103 au (    -57457.8115 kcal/mol)
# Free Energy Correction (Harmonic) : ZPVE + [H-TS]_T,R,V       =    -0.00944622 au (        -5.9276 kcal/mol)
# Gibbs Free Energy (Harmonic)      : E0 + ZPVE + [H-TS]_T,R,V  =   -91.57429725 au (    -57463.7391 kcal/mol)
# 

3
Iteration 17 Energy -91.56485103 (Optimized Structure)
 C   -0.5144063176   -0.0000000047    0.1645795165
 N    0.0609018094   -0.0000000048    1.2419920540
 H   -1.3594042072   -0.0000000048    1.0194684398

-1248.379728 ##Negative eigenvalue characterizing the TS
-0.038003 -0.000000  0.084134
 0.071381  0.000000 -0.012123
-0.539092  0.000000 -0.834035

 2104.774794
 0.362198 -0.000000  0.624794
-0.313716  0.000000 -0.554828
 0.043611 -0.000000  0.265124

 3070.885390
-0.051153 -0.000000  0.036968
-0.018025 -0.000000  0.004718
 0.859989  0.000000 -0.506051
```
