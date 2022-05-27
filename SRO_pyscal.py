#!/usr/bin/env python3
'''
################################################################################
# USAGE:: Run this script in a directory where SQS simulations are performed
# 
# USAGE:: pyscal to compute SRO parameters using fontain et al.
# https://pyscal.org/en/latest/_modules/pyscal/core.html#System.calculate_pmsro
################################################################################
'''

import sys, os
import pyscal as pc
import pyscal.crystal_structures as pcs
from ase import Atoms
from ase.io.vasp import write_vasp, read_vasp
import pandas as pd

gg={}
a=0; b=0; c=0;
for u in next(os.walk('.'))[1]:	
	print(u, end=' ')
	# Read initial and final states:
	hh = os.path.join(u, 'POSCAR_bestsqs')
	initial = read_vasp(hh)
	pos = initial.get_positions()
	ucell = initial.get_cell()
	
	type = initial.get_chemical_symbols()
	atom_type = sorted(list(set(type)))
	natm=len(atom_type)
	print ("Order of atom types:", atom_type)
	
	top = pc.System()
	top.read_inputfile(hh,format='poscar')
	#top.find_neighbors(method='cutoff', cutoff=1.2)
	top.find_neighbors(method='number', nmax=26)	
	#print("{:6s}  {:10.10s} {:10.10s} {:10.10s}".format( "Pairs", "1st_Shell", "2nd_Shell", "1st+2nd"  ) )
	for t in range(natm):
		for j in range(natm):
			sro = top.calculate_pmsro(reference_type=t+1, compare_type=j+1, average=True, shells=2, delta=True)
			if (j > t):
				#print("{}-{} {:9.5f} {:9.5f} {:9.5f}".format( atom_type[t], atom_type[j], sro[0], sro[1], sro[0]+sro[1]  ) )
				a=a+sro[0]
				b=b+sro[1]
				c=c+sro[0]+sro[1]
				gg[atom_type[t]+'-'+atom_type[j]] = [ round(sro[0],4), round(sro[1],4), sro[0]+sro[1] ]
			if (j == t):
				#print("{}-{} {:9.5f} {:9.5f} {:9.5f}".format( atom_type[t], atom_type[j], sro[0], sro[1], sro[0]+sro[1]   ) )
				a=a+sro[0]
				b=b+sro[1] 			
				c=c+sro[0]+sro[1]
				gg[atom_type[t]+'-'+atom_type[j]] = [ round(sro[0],4), round(sro[1],4), sro[0]+sro[1] ]
	df = pd.DataFrame(gg, index=['1st','2nd','1st+2nd']).T		
	print(df)
	df.to_csv(u+'.csv')
	print("{:^s}".format("-"*35))
	print("Tot_SRO of the system {:10.5f}".format(c))

