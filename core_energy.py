#!/usr/bin/env python3
	
'''
############################################################################
# USAGE  :: python3 sys.argv[0] 
# Author :: Asif Iqbal
# DATED  :: 01/12/2020
# NB     :: POSCAR can be in Cartesian/Direct coordinates.
# Calculates the core energy of the screw dislocations
# by sampling the different local chemical environment within a supercell
# by translating the atoms while keeping the screw dislocation fixed.
############################################################################
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
	pos = []; kk = []; lattice = []; sum = 0; dict = {}; g = 0
	file = open('POSCAR_perfect','r')
	firstline  = file.readline() # IGNORE first line comment
	alat = float( file.readline() )# scale
	Latvec1 = file.readline().split(); #print("{:9.6f} {:9.6f} {:9.6f}".format(float(Latvec1[0]),float(Latvec1[1]),float(Latvec1[2])))
	Latvec2 = file.readline().split(); #print("{:9.6f} {:9.6f} {:9.6f}".format(float(Latvec2[0]),float(Latvec2[1]),float(Latvec2[2])))
	Latvec3 = file.readline().split(); #print("{:9.6f} {:9.6f} {:9.6f}".format(float(Latvec3[0]),float(Latvec3[1]),float(Latvec3[2]))) 
	elementtype= file.readline(); #print ("{}".format(elementtype.split() ))
	atomtypes  = file.readline(); #print ("{}".format(atomtypes.split() ))
	Coordtype  = file.readline().split()
	if (Coordtype[0] == 'Direct' or Coordtype[0] == 'direct'): exit("First, Convert to Cartesian!")
	nat = [int(i) for i in atomtypes.split()]
	for i in nat: sum = sum + i; n_atoms = sum
	# print ("Number of atoms:", (n_atoms), end = '\n')	
	# Reading the Atomic positions				
	for x in range(int(n_atoms)):
		coord = [ float(i) for i in file.readline().split() ]
		pos = pos + [coord]
	file.close()
	#                              !!! This code index the atoms !!!
	#for index, line in enumerate(pos):
	#	dict[index] = line
	#for num, atm_num in dict.items():
	#	print("{} {}".format(num, atm_num) )
		#                      !!! TURN THIS ON IF YOU WANT A FILE IN XYZ FORMAT !!!
	#for j in range( len( elementtype.split() )):
	#	dict[elementtype.split()[j]] =  atomtypes.split()[j]; 
	#for l in dict:
	#	for k in range( int(dict[l]) ):
	#		#print( elementtype.split()[ math.floor( k/int(atomtypes.split()[0] ) ) ], pos[k][:] )
	#		print( l, pos[g][:] )
	#		g +=1
	return n_atoms,pos,firstline,alat,Latvec1,Latvec2,Latvec3,elementtype,atomtypes,Coordtype

def write_result(i,j,cnt,firstline,alat,Latvec1,Latvec2,Latvec3,elementtype,atomtypes,pos,n_atoms):
	with open( "POSTMP_"+str(cnt).zfill(2), 'w') as fdat1:
		fdat1.write( "{}\n".format( "POSTMP_"+str(i)+str(j) ) ) # Comment line in POSCAR
		fdat1.write( "{:5f}\n".format(alat) )
		fdat1.write( "{:15.11f} {:15.11f} {:15.11f}\n".format(float(Latvec1[0]),float(Latvec1[1]),float(Latvec1[2])) )
		fdat1.write( "{:15.11f} {:15.11f} {:15.11f}\n".format(float(Latvec2[0]),float(Latvec2[1]),float(Latvec2[2])) )
		fdat1.write( "{:15.11f} {:15.11f} {:15.11f}\n".format(float(Latvec3[0]),float(Latvec3[1]),float(Latvec3[2])) )
		fdat1.write( "{:5s}".format(elementtype) )
		fdat1.write( "{:5s}".format(atomtypes) )
		fdat1.write( "{:5s}\n".format(Coordtype[0]) )

 		# Displace the cell in the "X" direction.		
		for x in range(0, int(n_atoms), 1): 
			fdat1.write( "{:15.12f} {:15.12f} {:15.12f}\n".format(pos[x][0]+i,pos[x][1]+j,pos[x][2] ) )
	
# -------------------------------------- MAIN PROGRAM -------------------------------------- 
if __name__ == "__main__":
	n_atoms,pos,firstline,alat,Latvec1,Latvec2,Latvec3,elementtype,atomtypes,Coordtype = read_poscar();
	Ax = float(Latvec1[0]); 
	Ay = np.linalg.norm(Latvec2); 
	cnt = 0; cx = 0; cy = 0
	print("SYSTEM detected={}, #_atoms={}".format(elementtype.split(), n_atoms), end='\t\n' )
	print("Displaced in the X=<112> direction.")
	
	for i in np.linspace(0, Ax, 4, Ax):
		cx +=1
		for j in np.linspace(0, Ay, 4, Ay):
			cy += 1
			write_result(i,j,cnt,firstline,alat,Latvec1,Latvec2,Latvec3,elementtype,atomtypes,pos,n_atoms)
			L = 'dis_'+'X'+str(cx)+"_"+'Y'+str(cy)+'_'+str(cnt).zfill(2)
			shutil.rmtree( L , ignore_errors=True) #overwrite a directory
			
			os.mkdir( L )
			# copy the files to the directory.
			subprocess.call(['cp','-r','POSTMP_'+str(cnt).zfill(2), L ], shell = False)
			subprocess.call(['cp','-r','INCAR', L ], shell = False)
			subprocess.call(['cp','-r','POTCAR', L ], shell = False)
			subprocess.call(['cp','-r','KPOINTS', L ], shell = False)
			subprocess.call(['cp','-r','job.sh', L ], shell = False)
			subprocess.call(['cp','-r','input_dislo.babel', L ], shell = False)
			
			# Enter the directory.		
			os.chdir( L )
			subprocess.call(['cp','-r','POSTMP_'+str(cnt).zfill(2), 'CONTCAR'], shell = False)
			subprocess.call(['dislo', 'input_dislo.babel'], shell = False)
			subprocess.call(['cp','-r', 'POSCAR','../POSTMP_'+str(cnt).zfill(2)], shell = False)		
			os.chdir('../')		
			
			cnt += 1
			print ( "Ax,Ay = {:12.6f} {:12.6f} {:3d}".format( i,j, cnt ), end="\n" )
		
	print("DISPLACED FILES HAS BEEN GENERATED in the X=<112>, Y=<110> ... ")


			
	
	
