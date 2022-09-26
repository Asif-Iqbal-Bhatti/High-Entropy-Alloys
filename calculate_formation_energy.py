#!/usr/bin/env python3
'''
# AUTHOR:: ASIF IQBAL BHATTI
# USAGE:: TO COMPUTE THE FORMATION ENERGY OF STABLE alloys
# 
'''

import os, re
import pandas as pd
from ase.io import read
from ase.formula import Formula

inp_FILE1 = 'Ehull_lte_0.0_ARG2M'
inp_FILE2 = 'bulk_elementals'

out_File = inp_FILE1

subscript = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
elements_data, compos_data = {}, {}
for_ene = {}

for u in (next(os.walk(inp_FILE2))[1]):	
	vr = os.path.join(inp_FILE2, u, 'vasprun.xml')
	ene = read(vr, index=-1)
	t_ENE = ene.get_total_energy()
	t_ANU = ene.get_global_number_of_atoms()
	elements_data[u] = [t_ENE] + [t_ANU] + [t_ENE/t_ANU]

for u in (next(os.walk(inp_FILE1))[1]):	
	vr = os.path.join(inp_FILE1, u, '03_ISIF3snsb/vasprun.xml')
	ene = read(vr, index=-1)
	t_ENE = ene.get_total_energy()
	t_ANU = ene.get_global_number_of_atoms()
	t_FOR = Formula(u.split('_')[1]).count()
	compos_data[u.split('_')[1]] = [t_ENE] + [t_ANU] + [t_ENE/t_ANU]

#============ MAIN
for k, v in compos_data.items():
	a0 = re.findall('[A-Z][a-z]?', str(k))
	a = 0
	for i in a0:
		a += elements_data[i][2]
	for_ene[k] = v[2]-a
	
df = pd.DataFrame.from_dict(for_ene, orient='index')
df.to_csv(out_File+"-DATA.csv", header=['Ef/at'], index=True)
	
print(df)
