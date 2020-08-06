#!/usr/bin/env python3
#------------------------------------------------------------------------
# USAGE  :: python3 monte_carlo_HEA.py 
# Author :: Asif Iqbal
# DATED  :: 05/08/2020
# Metropolis Monte Carlo in a NVT (canonical) ensemble
# ADAPTED FROM:: https://chryswoods.com/intro_to_mc/part1/metropolis.html
# This script calculates the energy of the system by swapping the atoms
# and invoking the MC code to find the lowest structure.
#------------------------------------------------------------------------

import numpy as np
import os, sys, random, subprocess, shutil
import os.path, time

k = 8.617333262145E-5 # Boltzmann constant
T = 1000 # Temperature in Kelvin
sample = 10 # Number of sample could be # of atoms to swap

if os.path.exists('profile.dat'):
	os.remove('profile.dat') #this deletes the file
				
def read_poscar():
	pos = []; kk = []; lattice = []; sum = 0
	file = open('POSCAR','r')
	firstline  = file.readline() # IGNORE the comment line
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
	# print ("Number of atoms:", (n_atoms), end = '\n')	
	# Reading the Atomic positions				
	for x in range(int(n_atoms)):
		coord = file.readline().split()
		coord = [float(i) for i in coord]
		pos = pos + [coord]
	pos = np.array(pos)
	file.close()
	return n_atoms,pos,firstline,alat,Latvec1,Latvec2,Latvec3,elementtype,atomtypes,Coordtype

def calculate_energy():
	#old_energy = os.popen(" grep 'free  energy   TOTEN  =' OUTCAR | tail -1 | awk '{print $5 }' " ).read()
	#old_energy = float ( old_energy )
	f = open('OUTCAR',"r")
	lines = f.readlines()
	f.close()
	for i in lines:
		word = i.split()
		if "free  energy   TOTEN  =" in i:
			ii=lines.index(i)
	old_energy =  float (lines[ii].split()[4])
	return old_energy

def metropolis_MC(new_energy, old_energy, old_pos, new_pos, naccept, nreject):	
	tot_energy = []
	accept = False;
	# Accept if the energy goes down
	if (new_energy <= old_energy):
		accept = True
	else:
		# Apply the Monte Carlo test and compare
		# exp( -(E_new - E_old) / (k*T) ) >= rand(0,1)
		x = np.exp( -(new_energy - old_energy) / (k*T) )
		#print (x)
		if (x >= random.uniform(0.0,1.0)):
			accept = True
		else:
			accept = False
	if accept:
		# Accept the move
		naccept += 1; 
		print ("{} : {:10.6f}".format("Accept ratio", naccept/sample)  )
		tot_energy = new_energy
	else:
		# reject the move - restore the old coordinates
		nreject += 1
		print ("{} : {:10.6f}".format ( "Reject ratio", nreject/sample)  )
		new_pos[0] = old_pos[0]
		new_pos[1] = old_pos[1]
		new_pos[2] = old_pos[2]
		tot_energy = old_energy
		
	return new_pos, tot_energy, naccept, nreject

#------------------------------------MAIN PROGRAM--------------------------
#------------------------------------MAIN PROGRAM--------------------------
	
# First calculate the relaxed energy of the current SQS or SRO structure
naccept = 0; nreject = 0; 
old_energy = calculate_energy(); 
n_atoms, pos, firstline, alat, Latvec1,Latvec2,Latvec3, elementtype, atomtypes, Coordtype = read_poscar();
print ("{:20.30s} {:15.8f}".format('--> Initial system Energy', old_energy) )

for i in range(1, sample):
	old_pos = pos
	# Pick a random integer atom 1 (random.randint(a,b)) inclusive. Save the old coordinates
	rnd_atm1   = random.randint(0, n_atoms-1);
	#old_coords = ( pos[rnd_atm1][0], pos[rnd_atm1][1], pos[rnd_atm1][2] )
	#print ("Randomly selected atom", (old_coords), rnd_atm1-1, end = '\n')
	tmp_1 = pos[rnd_atm1][0]
	tmp_2 = pos[rnd_atm1][1]
	tmp_3 = pos[rnd_atm1][2]
	
	# Pick a random integer atom 2 (random.randint(a,b)) inclusive. save the old coordinates 
	rnd_atm2      = random.randint(0, n_atoms-1);
	#random_coords = ( pos[rnd_atm2][0], pos[rnd_atm2][1], pos[rnd_atm2][2] )
	#print ("Randomly selected atom", (random_coords), rnd_atm2-1, end = '\n')
	
	# Swapping the position of the two atoms and writing to a POSCAR file 

	pos[rnd_atm1][0] = pos[rnd_atm2][0]
	pos[rnd_atm1][1] = pos[rnd_atm2][1]
	pos[rnd_atm1][2] = pos[rnd_atm2][2]	
	pos[rnd_atm2][0] = tmp_1
	pos[rnd_atm2][1] = tmp_2
	pos[rnd_atm2][2] = tmp_3
	
	new_pos = pos
	print ("Swapping atoms", rnd_atm1-1, "and", rnd_atm2-1, end = '\t')
	
	#>>>>>>>>>-------- Writing a POSCAR File for new swapping -----------")
	fdata2 = open('POSCAR_'+str(i).zfill(3),'w')
	fdata2.write(firstline)
	fdata2.write("{:12.6f}\n".format(alat) )
	fdata2.write("{:15.12f} {:15.12f} {:15.12f}\n".format( float(Latvec1[0]),float(Latvec1[1]),float(Latvec1[2])) )
	fdata2.write("{:15.12f} {:15.12f} {:15.12f}\n".format( float(Latvec2[0]),float(Latvec2[1]),float(Latvec2[2])) )
	fdata2.write("{:15.12f} {:15.12f} {:15.12f}\n".format( float(Latvec3[0]),float(Latvec3[1]),float(Latvec3[2])) )
	fdata2.write(elementtype)
	fdata2.write(atomtypes)
	fdata2.write("{}\n".format( Coordtype[0] ) )
	for x in range(int(n_atoms)):
		fdata2.write("{:20.16f} {:20.16f} {:20.16f}\n".format(new_pos[x][0], new_pos[x][1], new_pos[x][2]) )	
	fdata2.close()
	#>>>>>>>>>-------- Ending a POSCAR File for new swapping -----------")
	
	shutil.rmtree('POS_'+str(i).zfill(3), ignore_errors=True) #overwrite a directory
	os.mkdir( 'POS_'+str(i).zfill(3) )
	
	subprocess.call(['cp','-r','POSCAR_'+str(i).zfill(3),'POS_'+str(i).zfill(3)], shell = False)
	subprocess.call(['cp','-r','INCAR','POS_'+str(i).zfill(3)], shell = False)
	subprocess.call(['cp','-r','POTCAR','POS_'+str(i).zfill(3)], shell = False)
	subprocess.call(['cp','-r','KPOINTS','POS_'+str(i).zfill(3)], shell = False)
	subprocess.call(['cp','-r','job.sh','POS_'+str(i).zfill(3)], shell = False)
	subprocess.call(['cp','-r','OUTCAR','POS_'+str(i).zfill(3)], shell = False)
	
	os.chdir('POS_'+str(i).zfill(3))
	shutil.copyfile('POSCAR_'+str(i).zfill(3), 'POSCAR' )
	
	subprocess.call(['sbatch','job.sh'], shell = False)	
	time.sleep(50)

	# Calculate the new energy of the swap atoms
	new_energy = calculate_energy();
	print ("Current system Energy:", (new_energy), end = '\t')
	
	new_pos, tot_energy, naccept, nreject = metropolis_MC(new_energy, old_energy, old_pos, new_pos, naccept, nreject)
	#old_energy = tot_energy
	#pos = new_pos

	os.chdir('../')
	
	with open('profile.dat', 'a') as fdata3:
		fdata3.write (" {:2d} {:15.8f} {:20.25s} {:12.8f} {:12.8f}\n".format(i, new_energy, 'POSCAR_'+str(i).zfill(3), naccept/sample, nreject/sample ))


	
	
