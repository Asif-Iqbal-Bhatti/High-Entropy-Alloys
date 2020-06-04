#!/usr/bin/env python3
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# USAGE :: ./python3 
# Metropolis Monte Carlo, in a NVT (canonical) ensemble, the Monte Carlo code is;
# https://chryswoods.com/intro_to_mc/part1/metropolis.html
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import numpy as np
import os, sys, random, subprocess, shutil
import 	os.path
import 	multiprocessing as mp

k = 8.617333262145E-5
T = 300
sample = 5

def read_poscar():
	pos = []; kk = []; lattice = []; sum = 0
	file = open('POSCAR','r')
	firstline  = file.readline() # IGNORE first line comment
	alat = float( file.readline() )# scale
	Latvec1 = file.readline().split(); #print("{:9.6f} {:9.6f} {:9.6f}".format(float(Latvec1[0]),float(Latvec1[1]),float(Latvec1[2])))
	Latvec2 = file.readline().split(); #print("{:9.6f} {:9.6f} {:9.6f}".format(float(Latvec2[0]),float(Latvec2[1]),float(Latvec2[2])))
	Latvec3 = file.readline().split(); #print("{:9.6f} {:9.6f} {:9.6f}".format(float(Latvec3[0]),float(Latvec3[1]),float(Latvec3[2]))) 
	elementtype= file.readline(); #print ("{}".format(elementtype.split() ))
	atomtypes  = file.readline(); #print ("{}".format(atomtypes.split() ))
	Coordtype  = file.readline().split()
	nat = atomtypes.split()
	nat = [int(i) for i in nat]
	for i in nat: sum = sum + i
	n_atoms = sum
	#print ("Number of atoms:", (n_atoms), end = '\n')	
	#>>>>>>>>>---------------Atomic positions------------------")				
	for x in range(int(n_atoms)):
		coord = file.readline().split()
		coord = [float(i) for i in coord]
		pos = pos + [coord]
	pos = np.array(pos)
	file.close()
	return n_atoms, pos, firstline, alat, Latvec1,Latvec2,Latvec3, elementtype, atomtypes, Coordtype

def calculate_energy():
	old_energy = os.popen(" grep 'free  energy   TOTEN  =' OUTCAR | tail -1 | awk '{print $5 }' " ).read()
	old_energy = float ( old_energy )
	return old_energy

def metropolis_MC(new_energy, old_energy, old_coords):	
	naccept = 0; nreject = 0;
	accept = False;
	# Automatically accept if the energy goes down
	if (new_energy <= old_energy):
			accept = True
	else:
			# Now apply the Monte Carlo test - compare
			# exp( -(E_new - E_old) / kT ) >= rand(0,1)
			x = np.exp( -(new_energy - old_energy) / (k*T) )
			if (x >= random.uniform(0.0,1.0)):
					accept = True
			else:
					accept = False
	
	if accept:
			# accept the move
			naccept += 1; print ("Acceptance:", naccept/sample) 
			total_energy = new_energy
	else:
			# reject the move - restore the old coordinates
			nreject += 1
			pos[atom_old][0] = old_coords[0]
			pos[atom_old][1] = old_coords[1]
			pos[atom_old][2] = old_coords[2]
			total_energy = old_energy
	return pos
#------------------------------------MAIN PROGRAM--------------------------
	
# First calculate the energy of the current SQS or SRO structure
old_energy   = calculate_energy();
n_atoms, pos, firstline, alat, Latvec1,Latvec2,Latvec3, elementtype, atomtypes, Coordtype = read_poscar();
print ("Energy of the system:", (old_energy), end = '\n')

for i in range(1, sample):	
	# Pick a random atom (random.randint(a,b) 
	# picks a random integer between a and b, inclusive. save the old coordinates
	atom_old   = random.randint(0, n_atoms-1);
	old_coords = ( pos[atom_old][0], pos[atom_old][1], pos[atom_old][2] )
	print ("Randomly selected atom", (old_coords), atom_old-1, end = '\n')
	tmp_1 = pos[atom_old][0]
	tmp_2 = pos[atom_old][1]
	tmp_3 = pos[atom_old][2]
	
	# Pick a random atom (random.randint(a,b) 
	# picks a random integer between a and b, inclusive)
	atom_rand     = random.randint(0, n_atoms-1);
	random_coords = ( pos[atom_rand][0], pos[atom_rand][1], pos[atom_rand][2] )
	print ("Randomly selected atom", (random_coords), atom_rand-1, end = '\n')
	
	# Swapping the position of the atoms and writing 
	# a POSCAR file and then submitting 
	# to VASP to calculate the new energy
	pos[atom_old][0] = pos[atom_rand][0]
	pos[atom_old][1] = pos[atom_rand][1]
	pos[atom_old][2] = pos[atom_rand][2]
	
	pos[atom_rand][0] = tmp_1
	pos[atom_rand][1] = tmp_2
	pos[atom_rand][2] = tmp_3
	
	#>>>>>>>>>-------- Writing a POSCAR File for new swapping -----------")
	fdata2 = open('POSCAR_'+str(i).zfill(3),'w')
	fdata2.write(firstline)
	fdata2.write("{:12.6f}\n".format(alat) )
	fdata2.write("{:15.12f} {:15.12f} {:15.12f}\n".format( float(Latvec1[0]),float(Latvec1[1]),float(Latvec1[2])) )
	fdata2.write("{:15.12f} {:15.12f} {:15.12f}\n".format( float(Latvec2[0]),float(Latvec2[1]),float(Latvec2[2])) )
	fdata2.write("{:15.12f} {:15.12f} {:15.12f}\n".format( float(Latvec3[0]),float(Latvec3[1]),float(Latvec3[2])) )
	fdata2.write(elementtype)
	fdata2.write(atomtypes)
	fdata2.write("{}\n".format( Coordtype[0]) )
	for x in range(int(n_atoms)):
		fdata2.write("{:20.16f} {:20.16f} {:20.16f}\n".format(pos[x][0], pos[x][1], pos[x][2]) )	
	fdata2.close()
	#shutil.copyfile("POSCAR_new", "POSCAR_"+str(i).zfill(3))
	
	# calculate the new energy of the swap atoms
	new_energy = calculate_energy();

	pos = metropolis_MC(new_energy, old_energy, old_coords)

	
	
	
