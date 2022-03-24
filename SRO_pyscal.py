#!/usr/bin/env python3.6
#
# USING pyscal to compute SRO parameters using fontain et al.
# https://pyscal.org/en/latest/_modules/pyscal/core.html#System.calculate_pmsro
#
import sys, os
import pyscal.core as pc
from ase.io.vasp import write_vasp, read_vasp
from tqdm import tqdm

top = pc.System()
cwd = os.path.abspath(os.path.dirname(__file__)) 
with open("SRO_column.dat", "w") as f:
	for k in next(os.walk('.'))[1]:
		d = os.path.join(cwd, k)

		top.read_inputfile(os.path.join(d, 'POSCAR_bestsqs'), format='poscar')
		top.find_neighbors(method='cutoff', cutoff=1.2)
		
		initial = read_vasp(os.path.join(d, 'POSCAR_bestsqs'))
		at_type = initial.get_chemical_symbols()
		atom_type = sorted(list(set(at_type)))
		print(k, atom_type,top.get_concentration())
		
		a=0; b=0; c=0; pairs_name = []; sro_1st_shell = []; sro_2st_shell = []; sro_1st_2nd_shell = []
		f.write("{}\n".format(k))
		f.write("{:6s}  {:10.10s} {:10.10s} {:10.10s}\n".format("Pairs", "1st Shell", "2nd Shell", "1st+2nd"))
		
		for t in range(5):
			for j in range(5):
				sro = top.calculate_pmsro(reference_type=t+1, compare_type=j+1, average=True, shells=2, delta=True)		
				if (j > t):
					f.write("{}-{} {:9.5f} {:9.5f} {:9.5f}\n".format(atom_type[t], atom_type[j], sro[0], sro[1], sro[0]+sro[1]))
					pairs_name.append(atom_type[t]+"-"+atom_type[j])
					sro_1st_shell.append(round(sro[0],4))
					sro_2st_shell.append(round(sro[1],4))
					sro_1st_2nd_shell.append(round(sro[0]+sro[1],4))
					a=a+sro[0]
					b=b+sro[1]
					c=c+sro[0]+sro[1]
				if (j == t):
					f.write("{}-{} {:9.5f} {:9.5f} {:9.5f}\n".format(atom_type[t], atom_type[j], sro[0], sro[1], sro[0]+sro[1]))
					pairs_name.append(atom_type[t]+"-"+atom_type[j])
					sro_1st_shell.append(round(sro[0],4))
					sro_2st_shell.append(round(sro[1],4))
					sro_1st_2nd_shell.append(round(sro[0]+sro[1],4))	
					a=a+sro[0]
					b=b+sro[1] 			
					c=c+sro[0]+sro[1]
	
		f.write("{}\n".format("-"*35))
		f.write ("{:5s} {:9.5f} {:9.5f} {:9.5f}\n".format('T_SRO', a, b, a+b))
		f.write("{}\n".format(pairs_name))
		f.write("{}\n".format(sro_1st_shell))
		f.write("{}\n".format(sro_2st_shell))
		f.write("{}\n".format(sro_1st_2nd_shell))
		f.write("{}\n".format("-"*50))
