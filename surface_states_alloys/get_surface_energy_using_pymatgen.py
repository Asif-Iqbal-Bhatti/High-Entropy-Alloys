#!/usr/bin/env python3

#https://www.materialssquare.com/blog/10-adsorption-energy-and-surface-energy-obtained-through-slab-structure

import os, sys
import pandas as pd
from pymatgen.io.vasp import Vasprun, Outcar
from pymatgen.core.structure import Molecule, Structure

cwd = os.getcwd()

fe_bcc = Outcar('../1_Fe_bcc/OUTCAR')
pure_H = Outcar('../3_Fe_atom/OUTCAR')
fe_ene = fe_bcc.final_fr_energy
h_ene  = pure_H.final_fr_energy
print(fe_ene, h_ene)

surf = {}
directories = [d for d in os.listdir(cwd) if os.path.isdir(os.path.join(cwd, d))]
for d in directories:
    se = Structure.from_file(d+'/CONTCAR')
    vpr = Outcar(d+'/OUTCAR')
    a = se.lattice.abc[0]
    b = se.lattice.abc[1]
    natoms = se.num_sites
    tot_ene = vpr.final_fr_energy
    
    surf_ene = ( tot_ene - natoms * (fe_ene/2.0) ) / 2.0*a*b
    #adsp_ene = (se.get_total_energy() - natoms * (fe_ene) - h_ene)
    surf[d] = [surf_ene]

df = pd.DataFrame(surf, dtype="float64[pyarrow]").T
df.columns = columns=['Surface Energy [eV/A**2]']
print(df)