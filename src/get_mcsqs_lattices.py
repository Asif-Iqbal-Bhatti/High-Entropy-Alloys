#!/usr/bin/env python3

#==================================================
# USAGE  :: AUTOMATE THE FINDING OF MCSQS LATTICE FROM LARGE DATA from VASP
# AUTHOR :: AIBEM2R
# DATED  :: 12/06/2024
#==================================================

import os
import sys
import pandas as pd
from tqdm import tqdm
from ase.io import read

lat_data = {}

inp_FILE1 = '../VASPRun_x0.20_800/VASPRun_x0.20_800/BEST'
outFILE = os.path.basename(inp_FILE1.rstrip('/'))

divisor = float(sys.argv[1]) # Change divisor based on input argument

for entry in tqdm(next(os.walk(inp_FILE1))[1]):
    try:
        pos = read(os.path.join(inp_FILE1, entry, "vasprun.xml"), index=-1)
        
        lavg = pos.cell.cellpar()[:3] / divisor 
        atm = len(pos)
        ene = pos.get_potential_energy() / atm
        vol = pos.get_volume()
        
        a, b, c = lavg
        a_prime = (a + b + c) / 3.0
        
        lat_data[entry] = [atm, a, b, c, a_prime, ene]
    except Exception as e:
        print(f"Error processing {entry}: {e}")
        continue

dfa = pd.DataFrame.from_dict(lat_data, orient='index')
dfa.columns = ['#atoms', '|a|', '|b|', '|c|', 'a_avg[A]', 'E[eV/atom]']
df_sorted = dfa.sort_values(by='E[eV/atom]')

min_index, max_index = df_sorted['E[eV/atom]'].idxmin(), df_sorted['E[eV/atom]'].idxmax()
min_row, max_row = df_sorted.loc[min_index], df_sorted.loc[max_index]

LE, HE = min_row['E[eV/atom]'], max_row['E[eV/atom]']

print(df_sorted, '\n')

avg_a_avg = dfa['a_avg[A]'].mean()
output = [
    f'AVG_over_all_data_sets,{avg_a_avg:6.4f}',
    f'LE,{min_index},{LE:6.4f}',
    f'HE,{max_index},{HE:6.4f}'
]

for line in output:
    print(line)

dfa.to_csv(f"{outFILE}.csv", header=True, index=True)

with open(f"{outFILE}.csv", "a") as f:
    f.write('\n')
    for line in output:
        f.write(line + '\n')

        
        
