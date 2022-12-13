#!/usr/bin/env python3
'''
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# AUTHOR:: AIB_EM
# USAGE PLOT DOS USING SUMO LIB
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'''
import os
import matplotlib.pyplot as pylt
from sumo.cli.dosplot import dosplot
from ase.formula import Formula
from ase.io import read

abs_path = os.path.abspath(os.path.join(os.getcwd(), '.'))

for n, u in enumerate(next(os.walk(abs_path))[1]):
	if not u.startswith('.') and n > 0:
		print(u, n)
		vprun = os.path.join(u, '03_DOS','vasprun.xml')
		tit = read(vprun, index=-1)

		w = f'{Formula(tit.get_chemical_formula()):latex}'
		v = f'{Formula(u.split("_")[1]):latex}'
		
		ff = dosplot(vprun, fonts='DejaVu Sans', ylabel=None, \
		yscale=1, zero_energy=None, zero_line=True, \
		height=6.0, width=8.0, plt=pylt, num_columns=1, \
		legend_cutoff=4)
		
		ff.title(v)
		ff.savefig( u+'.pdf', dpi=200, bbox_inches='tight')
		ff.clf()
		ff.close()
		del ff
