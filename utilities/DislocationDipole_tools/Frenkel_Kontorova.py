#!/usr/bin/env python3
	
'''
#------------------------------------------------------------------------
# USAGE  :: python3 frank_kontorova.py 
# Author :: Asif Iqbal
# DATED  :: 28/10/2020
# The multi-string Frenkel–Kontorova (MSFK) model.
# Python script to pull the strings in the 1/2<111> direction.
# """NOTA BINA: COORDINATES SHOULD BE IN CARTESIAN COORDINATES not DIRECT."""
# First, visualize the strings # near the dislocation line using the OVITO code.
# https://www.tandfonline.com/doi/pdf/10.1080/14786430903049104?needAccess=true
# In Ovito program atom numbering start from zero.
#------------------------------------------------------------------------
'''

import numpy as np
import os, sys, random, subprocess, shutil
import matplotlib as mp
from ase.constraints import FixAtoms, FixedPlane, FixedLine
from ase import Atoms
import ase.io
from ase.io import write, read
from ase.io.vasp import write_vasp, read_vasp

def read_poscar():
	pos = []
	kk = []
	lattice = []
	sum = 0
	dict = {}
	g = 0
	with open('POSCAR','r') as file:
		firstline  = file.readline() # IGNORE first line comment
		alat = float( file.readline() )# scale
		Latvec1 = file.readline().split()
		Latvec2 = file.readline().split()
		Latvec3 = file.readline().split()
		elementtype= file.readline()
		atomtypes  = file.readline()
		Coordtype  = file.readline().split()
		if Coordtype[0] in ['Direct', 'direct']: exit("Convert to Cartesian!")
		nat = [int(i) for i in atomtypes.split()]
		for i in nat: sum = sum + i; n_atoms = sum
		for _ in range(int(n_atoms)):
			coord = [ float(i) for i in file.readline().split() ]
			pos = pos + [coord]
		pos = [ [float(pos[j][i]) for i in range(3)] for j in range(n_atoms) ]
	return n_atoms,pos,firstline,alat,Latvec1,Latvec2,Latvec3,elementtype,atomtypes,Coordtype

def write_result(d,cnt,firstline,alat,Latvec1,Latvec2,Latvec3,elementtype,atomtypes,pos,n_atoms):
	with open(f"POSCAR_{str(cnt).zfill(2)}", 'w') as fdat1:
		fdat1.write(f"POSCAR_{str(d)}\n")
		fdat1.write( "{:5f}\n".format(alat) )
		fdat1.write( "{:15.11f} {:15.11f} {:15.11f}\n".format(float(Latvec1[0]),float(Latvec1[1]),float(Latvec1[2])) )
		fdat1.write( "{:15.11f} {:15.11f} {:15.11f}\n".format(float(Latvec2[0]),float(Latvec2[1]),float(Latvec2[2])) )
		fdat1.write( "{:15.11f} {:15.11f} {:15.11f}\n".format(float(Latvec3[0]),float(Latvec3[1]),float(Latvec3[2])) )
		fdat1.write( "{:5s}".format(elementtype) )
		fdat1.write( "{:5s}".format(atomtypes) )
		fdat1.write( "{:5s}\n".format(Coordtype[0]) )
# "x" is a string number. Choose the strings # (Atom indexes) from the OVITO program.
		for x in range(int(n_atoms)): 
      #                            ''' !!! 7 strings => 14 atoms !!! '''
			if x not in [
				160,
				161,
				170,
				259,
				346,
				354,
				438,
				293,
				294,
				305,
				397,
				212,
				386,
				213,
			]: 
				fdat1.write( "{:12.9f} {:12.9f} {:12.9f}\n".format(pos[x][0],pos[x][1],pos[x][2] ) )
			elif (x == 160):
				fdat1.write("{:12.9f} {:12.9f} {:12.9f}\n".format(pos[x][0],pos[x][1],pos[x][2]+d ) )
			elif (x == 161):
				fdat1.write("{:12.9f} {:12.9f} {:12.9f}\n".format(pos[x][0],pos[x][1],pos[x][2]+d ) )
			elif (x == 170):
				fdat1.write("{:12.9f} {:12.9f} {:12.9f}\n".format(pos[x][0],pos[x][1],pos[x][2]+d ) )
			elif (x == 259):
				fdat1.write("{:12.9f} {:12.9f} {:12.9f}\n".format(pos[x][0],pos[x][1],pos[x][2]+d ) )
			elif (x == 346):
				fdat1.write("{:12.9f} {:12.9f} {:12.9f}\n".format(pos[x][0],pos[x][1],pos[x][2]+d ) )
			elif (x == 354):
				fdat1.write("{:12.9f} {:12.9f} {:12.9f}\n".format(pos[x][0],pos[x][1],pos[x][2]+d ) )
			elif (x == 438):
				fdat1.write("{:12.9f} {:12.9f} {:12.9f}\n".format(pos[x][0],pos[x][1],pos[x][2]+d ) )
			elif (x == 293):
				fdat1.write("{:12.9f} {:12.9f} {:12.9f}\n".format(pos[x][0],pos[x][1],pos[x][2]+d ) )
			elif (x == 294):
				fdat1.write("{:12.9f} {:12.9f} {:12.9f}\n".format(pos[x][0],pos[x][1],pos[x][2]+d ) )
			elif (x == 305):
				fdat1.write("{:12.9f} {:12.9f} {:12.9f}\n".format(pos[x][0],pos[x][1],pos[x][2]+d ) )
			elif (x == 397):
				fdat1.write("{:12.9f} {:12.9f} {:12.9f}\n".format(pos[x][0],pos[x][1],pos[x][2]+d ) )
			elif (x == 212):
				fdat1.write("{:12.9f} {:12.9f} {:12.9f}\n".format(pos[x][0],pos[x][1],pos[x][2]+d ) )
			elif (x == 386):
				fdat1.write("{:12.9f} {:12.9f} {:12.9f}\n".format(pos[x][0],pos[x][1],pos[x][2]+d ) )
			else:
				fdat1.write("{:12.9f} {:12.9f} {:12.9f}\n".format(pos[x][0],pos[x][1],pos[x][2]+d ) )

# -------------------------------------- MAIN PROGRAM -------------------------------------- 

if __name__ == "__main__":
	n_atoms,pos,firstline,alat,Latvec1,Latvec2,Latvec3,elementtype,atomtypes,Coordtype = read_poscar();
	Burgers = float(Latvec3[2])/2
	p_constr = []

	print(f"SYSTEM detected={elementtype.split()}, #_atoms={n_atoms}", end='\t\n')
	print("Displaced by 1 Burgers vector along the <111>.")
	atom_Freeze = [ 160, 161, 170, 259, 346, 354, 438, 293, 294, 305, 397, 212, 386, 213 ]

	for cnt, d in enumerate(np.linspace(0, Burgers, 15, Burgers)):
		write_result(d,cnt,firstline,alat,Latvec1,Latvec2,Latvec3,elementtype,atomtypes,pos,n_atoms)
		shutil.rmtree(f'POS_{str(cnt).zfill(2)}', ignore_errors=True)
		os.mkdir(f'POS_{str(cnt).zfill(2)}')

		subprocess.call(
			['cp', '-r', f'POSCAR_{str(cnt).zfill(2)}', f'POS_{str(cnt).zfill(2)}'],
			shell=False,
		)

		subprocess.call(
			['cp', '-r', 'INCAR', f'POS_{str(cnt).zfill(2)}'], shell=False
		)

		subprocess.call(
			['cp', '-r', 'POTCAR', f'POS_{str(cnt).zfill(2)}'], shell=False
		)

		subprocess.call(
			['cp', '-r', 'KPOINTS', f'POS_{str(cnt).zfill(2)}'], shell=False
		)

		subprocess.call(
			['cp', '-r', 'job.sh', f'POS_{str(cnt).zfill(2)}'], shell=False
		)


		os.chdir(f'POS_{str(cnt).zfill(2)}')

		#----------------- Constraining the atoms ---------------------
		#shutil.copyfile( 'POSCAR', 'POSCAR_'+str(cnt).zfill(2) )
		initial = read(f'POSCAR_{str(cnt).zfill(2)}')
		p_pos   = initial.get_positions()
		p_ucell = initial.get_cell()
		p_type  = initial.get_chemical_symbols()
		p_ll = Atoms(positions=p_pos, symbols=p_type, cell=p_ucell, pbc=(1,1,1) )

		### TO FIX ALL THE ATOMS
		#for atom in p_ll:
		#	p_constr.append( FixedPlane(atom.index, ( 0, 0, 1)) )
		#p_constr = [ FixedPlane( atom.index, ( 0, 0, 1) )  for atom in p_ll ]

		#### TO FIX ONLY CERTAIN ATOMS
		p_constr.extend(
			FixAtoms(indices=[atom.index for atom in p_ll if atom.index == item])
			for item in atom_Freeze
		)

		p_ll.set_constraint(p_constr)

		write_vasp("POSCAR", atoms=p_ll, direct=False, vasp5=True, ignore_constraints=False)

		#subprocess.call(['../selective.sh'], shell = False)
		subprocess.call(
			['cp', '-r', 'POSCAR', f'../POSCAR_{str(cnt).zfill(2)}'], shell=False
		)

		os.chdir('../')
		print ( "b = {:6.6f} {:3d}".format( d, cnt ), end="\n" )
	print("DISPLACED FILES HAS BEEN GENERATED ... ")


			
