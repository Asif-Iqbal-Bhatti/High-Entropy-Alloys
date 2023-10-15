#!/usr/bin/env python3

'''
#***************************************************************************************
# File: sys.argv[0]
# Author: Asifem2r/AIB_EM
# Date: 2023-10-13
# Description: To shake or rattle atom positions and nudge lattice vectors
# 
# Copyright (c) 2023 Asif Iqbal
# 
# Shake a system to avoid the minima
# atoms.rattle(stdev=0.001, seed=9875468, rng=None)
# write('POSCAR_reordered', atoms, format='vasp', direct=False)
# 
#***************************************************************************************/
'''

import random
import numpy as np
from ase.io import read, write

atoms = read('CONTCAR')
symbols = atoms.get_chemical_symbols()

max_displacement_atoms = 0.01  # Adjust as needed
max_displacement_lattice = 1.0  # Adjust as needed

# Shake the atomic positions
for atom in atoms:
    displacement = np.random.uniform(-max_displacement_atoms, max_displacement_atoms, size=3)
    atom.position += displacement

# Shake the lattice parameters
cell = atoms.get_cell()
size=(3, 3)
identity = np.eye(size[0])
displacement = np.random.uniform(-max_displacement_lattice, max_displacement_lattice, size=size)
displacement = displacement * identity
new_cell = cell + displacement
atoms.set_cell(new_cell)

# Write the modified structure to a new POSCAR file
output_poscar_file = 'POSCAR'
write(output_poscar_file, atoms, format='vasp', direct=True)
print(f"Modified structure written to {output_poscar_file}")
