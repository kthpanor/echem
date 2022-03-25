# Running a simulation

Workflow with example in openMM.

Can also link to good tutorial with [lysozyme](http://www.mdtutorials.com/gmx/lysozyme/index.html) etc (gromacs).

```python

```

## AIMD example with XTB

See https://xtb-docs.readthedocs.io/en/latest/md.html for more details.

- Want to run from notebook
    - Run without creating all files

- Visualize trajectory

```python
import numpy as np
import matplotlib.pyplot as plt

def write_xyz(mol_str,fname):
    input_str = '$coord\n'
    for i in np.arange(len(mol_str.split('\n'))):
        if mol_str.split('\n')[i]:
            data = mol_str.split('\n')[i].split(' ')
            tmp_str = ''
            for i in np.arange(len(data)):
                if data[i]:
                    if data[i].isalpha():
                        save_el = data[i].lower()
                    else:
                        tmp_str += data[i] + '   '
            tmp_str += save_el
            input_str += tmp_str + '\n'
    input_str += '$end'
    f = open(fname, 'w'); f.write(input_str); f.close()
    return False
```

```python
md_str = '''$md
   temp=300
   time= 25.0
   dump=100.0
   step=  2.5
$end
'''

mol_str = '''
O       0.0000000000     0.0000000000     0.1178336003
H      -0.7595754146    -0.0000000000    -0.4713344012
H       0.7595754146     0.0000000000    -0.4713344012
'''


mol_name = 'molecule'
md_name = 'md.inp'
save_xyz(mol_str,mol_name)
f = open('md.inp', 'w'); f.write(md_str); f.close()
```

```python
# !xtb $mol_name --input $md_name --omd
```

```python
traj = open('xtb.trj', 'r').read()
energies = []
xyzs = []

snapshots = traj.split('energy:')[1:]
n_atoms = traj.split('\n')[0]
for i in np.arange(len(snapshots)):
    energies.append(np.float((snapshots[i].split(' ')[1])))
    tmp_snapshot = snapshots[i].split(')\n')[1]
    xyzs.append(tmp_snapshot.split('\n{}\n'.format(n_atoms))[0])
    
plt.figure()
plt.plot(energies)
plt.show()
```
