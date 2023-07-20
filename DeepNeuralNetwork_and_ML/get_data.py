#!/usr/bin/env python3

import os, sys
import pandas as pd
import numpy as np
from ase.io import read, write
import dpdata


inpFilename = 'OUTCAR'
outFilename = 'heaData.extxyz'
os.system('rm -r ' + outFilename)

#
data = dpdata.LabeledSystem('OUTCAR', fmt = 'vasp/outcar') 
print(f'# the data contains {len(data)} frames, {data["coords"].shape}')

#
all_atoms = read(inpFilename, format='vasp-out', index=':')
print("READING LATTICE DIM:: ", all_atoms[-1].get_positions().shape)
print("# of Frames:: ", len(all_atoms))

for i, curr_atoms in enumerate(all_atoms): 
  write(outFilename, curr_atoms, append=True, format='extxyz')

# CHECK
dd = []
rd = read(outFilename, format='extxyz', index=':')
for k, i in enumerate(rd):
    dd.append({'Frame': k,
                'Formula':i.get_chemical_formula(),
                'TotENE':i.get_potential_energy()})

df = pd.DataFrame(dd)
print(df)
minENEinMD = df['TotENE'].idxmin()

print('Min Frame idx and Frame TotENE:: ', minENEinMD, df['TotENE'].min())
print(all_atoms[minENEinMD].get_potential_energy())

write('CONTCAR_min', all_atoms[minENEinMD], format='vasp')