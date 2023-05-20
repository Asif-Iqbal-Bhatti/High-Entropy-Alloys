#!/usr/bin/env python3

import numpy as np
from ase.io import read, write


inpFilename = 'OUTCAR'
outFilename = 'nequip-data.extxyz'

# read all frames into a list of ase.Atoms objects
all_atoms = read(inpFilename, format='vasp-out', index=':')
print(all_atoms[-1].get_positions().shape)

for i, curr_atoms in enumerate(all_atoms): 
  write(outFilename, curr_atoms, append=True, format='extxyz')