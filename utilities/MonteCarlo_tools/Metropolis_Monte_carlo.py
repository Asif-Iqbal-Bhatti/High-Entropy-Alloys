#!/usr/bin/env python3
'''
#------------------------------------------------------------------------
# USAGE  :: PYTHON3 SRO_vs_ENERGY.py
# AUTHOR :: @ASIFIQBAL
# DATED  :: 02/07/2022
# METROPOLIS MONTE CARLO IN A NVT (CANONICAL) ENSEMBLE
# ADAPTED FROM:: HTTPS://CHRYSWOODS.COM/INTRO_TO_MC/PART1/METROPOLIS.HTML
# THIS SCRIPT CALCULATES THE ENERGY OF THE SYSTEM BY SWAPPING THE ATOMS
# AND INVOKING THE MC CODE TO FIND THE LOWEST STRUCTURE.
#------------------------------------------------------------------------
'''

import numpy as np
import os, sys, random, subprocess, shutil
import os.path, time

k = 8.617333262145E-5 # Boltzmann constant
T = 500 # Temperature in Kelvin
sample = 200 # Number of sample could be # of atoms to swap

out_files = ["profile.csv", "accept.dat"];
for f in out_files: 
	if os.path.exists(f):
		os.remove(f) #this deletes the file
				
def read_poscar():
	pos = []
	kk = []
	lattice = []
	sum = 0
	with open('CONTCAR','r') as file:
		firstline  = file.readline() # IGNORE first line comment
		alat = float( file.readline() )# scale
		Latvec1 = file.readline().split()
		Latvec2 = file.readline().split()
		Latvec3 = file.readline().split()
		elementtype= file.readline()
		atomtypes  = file.readline()
		Coordtype  = file.readline().split()
		nat = atomtypes.split()
		nat = [int(i) for i in nat]
		for i in nat: sum = sum + i
		n_atoms = sum
			# print ("Number of atoms:", (n_atoms), end = '\n')	
			# Reading the Atomic positions				
		for _ in range(int(n_atoms)):
			coord = file.readline().split()
			coord = [float(i) for i in coord]
			pos = pos + [coord]
		pos = np.array(pos)
	return n_atoms,pos,firstline,alat,Latvec1,Latvec2,Latvec3,elementtype,atomtypes,Coordtype

def calculate_energy():
	with open('OUTCAR',"r") as f:
		lines = f.readlines()
	for i in lines:
		word = i.split()
		if "free  energy   TOTEN  =" in i:
			ii=lines.index(i)
	return float(lines[ii].split()[4])

def metropolis_MC(new_energy, old_energy, naccept, nreject):	
	a_energy = []
	r_energy = []
	accept = False;
	# Accept if the energy goes down
	if (new_energy <= old_energy):
		accept = True
	else:
		# Apply the Monte Carlo test and compare
		# exp( -(E_new - E_old) / kT ) >= rand(0,1)
		x = np.exp( -(new_energy - old_energy) / (k*T) )
		#print (x)
		accept = x >= random.uniform(0.0,1.0)
	if accept:
		# Accept the move
		naccept += 1; 
		#print ("{}: {:10.3f}%".format("Accept ratio", (naccept/sample)*100  )  )
		a_energy = new_energy
		yes="Accept"
	else:
		# reject the move - restore the old coordinates
		nreject += 1
		#print ("{}: {:10.3f}%".format ("Reject ratio", (nreject/sample)*100 )  )
		r_energy = old_energy	
		yes="Reject"
	return a_energy, r_energy, naccept, nreject, yes

def write_result(i,new_energy,old_energy,SRO,naccept,nreject,sample,yes):
	with open('profile.csv', 'a') as fdata3:
		fdata3.write(
			"{:3d}, {:9.10s}, {:11.8f}, {:8.5f}, {:8.3f}, {:8.3f}, {:s}\n".format(
				i,
				f'POS_{str(i).zfill(3)}',
				new_energy,
				SRO,
				(naccept / sample) * 100,
				(nreject / sample) * 100,
				yes,
			)
		)


	if (yes=='Accept'):
		with open('accept.dat', 'a') as fdata4:
			fdata4.write(
				"{:9.10s}, {:11.6f}, {:8.5f}, {:8.3f}, {:8.3f}, {:s}\n".format(
					f'POS_{str(i).zfill(3)}',
					new_energy,
					SRO,
					(naccept / sample) * 100,
					(nreject / sample) * 100,
					yes,
				)
			)
			
#=====================  MAIN PROGRAM
if __name__ == "__main__":
	# FIRST OBTAIN THE GROUND/OPTIMIZED ENERGY OF 
	# THE CURRENT SQS IN AN IDEAL POSITION
	naccept = 0
	nreject = 0;
	old_energy = calculate_energy();
	n_atoms, pos, firstline, alat, Latvec1,Latvec2,Latvec3, elementtype, atomtypes, Coordtype = read_poscar();
	print("Starting simulation energy: {:15.8f}".format(old_energy), end = '\n')
	print("T={:4f} K Sample={:5d} Atoms={:5d}\n".format(T, sample, n_atoms))

	with open('profile.csv', 'a') as fdata3:
		fdata3.write ("{:17s}, {:17s}, {:11.12s}, {:6.6s}, {:6s}, {:8s}\n".format(" ", " ","Ediff", "SRO", "Accept[%]", "Reject[%]" ))

	SRO=1.0
	for i in range(1, sample):

		os.chdir(f'POS_{str(i).zfill(3)}')
		new_energy = calculate_energy()
		#print('{:3d} Energy in POS_{:3s}: {:15.6f} {:13.6f}'.format(i, str(i).zfill(3), new_energy, SRO), end = '\t')
		a_energy, r_energy, naccept, nreject, yes = metropolis_MC(new_energy, old_energy, naccept, nreject)
		#print (old_energy, new_energy)

		os.chdir('../')

		write_result(i,new_energy,old_energy,SRO,naccept,nreject,sample,yes)
		old_energy = new_energy

	#print('Accepted:: {:3d}, Rejected:: {:3d}'.format(naccept, nreject), end = '\n')
	with open('profile.csv', 'a') as fdata3:
		fdata3.write ('Accepted:: {:3d}, Rejected:: {:3d}\n'.format(naccept, nreject) )
		
	
