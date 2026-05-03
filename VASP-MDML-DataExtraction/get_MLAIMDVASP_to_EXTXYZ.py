#!/usr/bin/env python3

# AUTHOR:: AsifEM2R

import json, os, sys
from tqdm import tqdm
from ase import Atoms
import numpy as np
from py4vasp import Calculation
from monty.json import MontyEncoder
from ase.io import write
from ase.calculators.singlepoint import SinglePointCalculator


# --- PATH SETUP ---
EXTXYZ_NAME = 'mlmd_from_outcar.extxyz'
xyzpath = os.path.join(EXTXYZ_NAME)
if os.path.exists(xyzpath):
    os.remove(xyzpath)

    
calc = Calculation.from_file("vaspout.h5")
print(dir(calc.energy[:]))
print(dir(calc))

def MLMD_AIMD_to_EXTXYZ(xyzpath, calc, buffer_size=100):    
    #struct_tomdtraj = calc.structure[:].to_mdtraj()
    #print(dir(struct_tomdtraj))
    #struct_tomdtraj.save_xyz('mdtraj_data.xyz')
    
    md_steps = calc.structure[:].number_steps()
    name = calc.topology
    buffer = []

    for i in tqdm(range(md_steps), desc="Processing MLAIMD frames"):
        atoms = calc.structure[i].to_ase()
        energy = calc.energy[i].read()['total energy   ETOTAL']
        forces = calc.force[i].read()['forces']
        stress = calc.stress[i].read()['stress']

        stress_ASE = np.array(
            [
                stress[0][0],
                stress[1][1],
                stress[2][2],
                stress[1][2],
                stress[2][0],
                stress[0][1],
            ]
        )
        # internal stress
        stress_eVA3 = -1 * stress_ASE / 1602.1766208  # from kbar to eV/Angstrom^3
                    
        atoms.calc = SinglePointCalculator(atoms, energy=energy, forces=forces, stress=stress_eVA3)
        atoms.info = {"System": f"{name}"}
        buffer.append(atoms)

        if len(buffer) >= buffer_size:
            write(f'{xyzpath}.extxyz', buffer, append=True)
            buffer.clear()

    if buffer:
        write(f'{xyzpath}', buffer, append=True)

######################################################

MLMD_AIMD_to_EXTXYZ(EXTXYZ_NAME, calc)
