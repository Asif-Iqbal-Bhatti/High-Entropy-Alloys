#!/usr/bin/env python3.6
#
# USING pyscal to compute SRO parameters using fontain et al.
# https://pyscal.org/en/latest/_modules/pyscal/core.html#System.calculate_pmsro
#
import sys
import pyscal as pc
import pyscal.crystal_structures as pcs
from ase import Atoms
from ase.io.vasp import write_vasp, read_vasp

# Read initial and final states:
initial = read_vasp(sys.argv[1])
pos = initial.get_positions()
ucell = initial.get_cell()

type = initial.get_chemical_symbols()
atom_type = sorted(list(set(type)))
print ("Order of the atom types", atom_type)

top = pc.System()
top.read_inputfile(sys.argv[1],format='poscar')
top.find_neighbors(method='cutoff', cutoff=1.2)

a=0; b=0; c = 0; pairs_name = []; sro_1st_shell = []; sro_2st_shell = []; sro_1st_2nd_shell = []
print("{:6s}  {:10.10s} {:10.10s} {:10.10s}".format( "Pairs", "1st Shell", "2nd Shell", "1st+2nd"  ) )
for t in range(5):
	for j in range(5):
		sro = top.calculate_pmsro(reference_type=t+1, compare_type=j+1, average=True, shells=2, delta=True)
		if (j > t):
			print("{}-{}, {:9.5f} {:9.5f} {:9.5f}".format( atom_type[t], atom_type[j], sro[0], sro[1], sro[0]+sro[1]  ) )
			pairs_name.append(atom_type[t]+"-"+atom_type[j])
			sro_1st_shell.append(sro[0])
			sro_2st_shell.append(sro[1])
			sro_1st_2nd_shell.append(sro[0]+sro[1])
			a=a+sro[0]
			b=b+sro[1]
			c=c+sro[0]+sro[1]
		if (j == t):
			print("{}-{}, {:9.5f} {:9.5f} {:9.5f}".format( atom_type[t], atom_type[j], sro[0], sro[1], sro[0]+sro[1]   ) )
			pairs_name.append(atom_type[t]+"-"+atom_type[j])
			sro_1st_shell.append(sro[0])
			sro_2st_shell.append(sro[1])
			sro_1st_2nd_shell.append(sro[0]+sro[1])			
			a=a+sro[0]
			b=b+sro[1] 			
			c=c+sro[0]+sro[1]

print ("Total SRO of the system {:9.5f} 1st+2nd shell".format(a+b))
print (pairs_name)
print (sro_1st_shell)
print (sro_2st_shell)
print (sro_1st_2nd_shell)
