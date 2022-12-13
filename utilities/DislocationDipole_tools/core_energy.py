#!/usr/bin/env python3
	
'''
############################################################################
# USAGE  :: python3 sys.argv[0] 
# Author :: Asif Iqbal
# DATED  :: 13/12/2020
# NB     :: POSCAR should be in Cartesian coordinates.
# Calculate the core energy of the screw dislocations
# by sampling different local chemical environment within a supercell
# by translating the atoms keeping the screw dislocation fixed.
############################################################################
'''

import numpy as np
import os, sys, random, subprocess, shutil

def read_POSCAR():
	pos = []
	kk = []
	lattice = []
	sum = 0;
	with open('POSCAR_perfect','r') as file:
		firstline  = file.readline() # IGNORE first line comment
		alat = float( file.readline() )# scale
		Latvec1 = file.readline().split()
		Latvec2 = file.readline().split()
		Latvec3 = file.readline().split()
		elementtype= file.readline()
		atomtypes  = file.readline()
		Coordtype  = file.readline().split()
		if Coordtype[0] in ['Direct', 'direct']: exit("First, Convert to Cartesian!")
		nat = [int(i) for i in atomtypes.split()]
		for i in nat: sum = sum + i; n_atoms = sum
		for _ in range(int(n_atoms)):
			coord = [ float(i) for i in file.readline().split() ]
			pos = pos + [coord]
	return n_atoms,pos,firstline,alat,Latvec1,Latvec2,Latvec3,elementtype,atomtypes,Coordtype

def write_result(i,j,cnt,firstline,alat,Latvec1,Latvec2,Latvec3,elementtype,atomtypes,pos,n_atoms):
	with open(f"iniTMP_{str(cnt).zfill(2)}", 'w') as fdat1:
		fdat1.write(f"iniTMP_{str(i)}_{str(j)}\n")
		fdat1.write( "{:5f}\n".format(alat) )
		fdat1.write( "{:15.11f} {:15.11f} {:15.11f}\n".format(float(Latvec1[0]),float(Latvec1[1]),float(Latvec1[2])) )
		fdat1.write( "{:15.11f} {:15.11f} {:15.11f}\n".format(float(Latvec2[0]),float(Latvec2[1]),float(Latvec2[2])) )
		fdat1.write( "{:15.11f} {:15.11f} {:15.11f}\n".format(float(Latvec3[0]),float(Latvec3[1]),float(Latvec3[2])) )
		fdat1.write( "{:5s}".format(elementtype) )
		fdat1.write( "{:5s}".format(atomtypes) )
		fdat1.write( "{:5s}\n".format(Coordtype[0]) )
 		# Displace the cell in the "X" & "Y" direction.
		for x in range(int(n_atoms)): 
			fdat1.write( "{:15.12f} {:15.12f} {:15.12f}\n".format(pos[x][0]+i,pos[x][1]+j,pos[x][2] ) )
	
# -------------------------------------- MAIN PROGRAM -------------------------------------- 
if __name__ == "__main__":
	n_atoms,pos,firstline,alat,Latvec1,Latvec2,Latvec3,elementtype,atomtypes,Coordtype = read_POSCAR();
	Ax = float(Latvec1[0]);
	Ay = np.linalg.norm(Latvec2[1]);
	cnt = 0
	cx = 0
	cy = 0
	b = 3.40/2 * np.linalg.norm([1,1,1])
	a = 2*b /np.sqrt(3)
	with open("FILE_INFO.yaml",'w') as file:
		for i in np.arange(0, Ax, np.sqrt(2/3)*3.40):
			for j in np.arange(0, Ay, np.sqrt(2/3)*3.40):
				write_result(i,j,cnt,firstline,alat,Latvec1,Latvec2,Latvec3,elementtype,atomtypes,pos,n_atoms)
				L = 'dis_'+'X'+str(cx).zfill(2)+"_"+'Y'+str(cy).zfill(2)+'_'+str(cnt).zfill(2)
				shutil.rmtree( L , ignore_errors=True ) #overwrite a directory

				os.mkdir( L )

							# copy the files to the directory.
				subprocess.call(['cp', '-r', f'iniTMP_{str(cnt).zfill(2)}', L], shell = False)
				subprocess.call(['cp','-r','INCAR', L ], shell = False)
				subprocess.call(['cp','-r','POTCAR', L ], shell = False)
				subprocess.call(['cp','-r','KPOINTS', L ], shell = False)
				subprocess.call(['cp','-r','input_dislo.babel', L ], shell = False)

				# Enter the directory.		
				os.chdir( L )

				with open("job.sh", 'w') as job:
					job.write("#!/bin/bash\n")
					job.write(f"#SBATCH -J core_{cnt}\n")
					job.write("#SBATCH -o slurm.out\n")
					job.write("#SBATCH -e slurm.err\n")
					job.write("#SBATCH -t 50:40:00\n")
					job.write("#SBATCH --nodes=6\n")
					job.write("#SBATCH --ntasks-per-node=40\n")
					job.write("#SBATCH --partition=COMPUTE-SHORT\n")
					job.write("module load intel/2018.4\n")
					job.write("module load intel/2018.4/openmpi/3.1.4\n")
					job.write("srun --resv-ports --mpi=pmi2 /softs/VASP544-INTEL2018-OMPI3.1.4-TBDYN/vasp.5.4.4/bin/vasp_std\n")
				subprocess.call(
					['cp', '-r', f'iniTMP_{str(cnt).zfill(2)}', 'CONTCAR'], shell=False
				)

				subprocess.call(['dislo', 'input_dislo.babel'], shell = False)
				K = 'disFIN_'+'X'+str(cx).zfill(2)+"_"+'Y'+str(cy).zfill(2)+'_'+str(cnt).zfill(2)
				subprocess.call( ['cp','-r', 'POSCAR', K ], shell = False)
				subprocess.call(
					['cp', '-r', K, f'../iniTMP_{str(cnt).zfill(2)}'], shell=False
				)

				os.chdir('../')		

				file.write ( "Ax, Ay = {:12.6f},{:12.6f} {:3d}\n".format( i,j, cnt ) )

				cnt += 1
				cy  += 1
			cy = 0	# reset counter
			cx +=1

		print("DISPLACED FILES HAS BEEN GENERATED in the X=<112>, Y=<110> ... ")
		print(f"SYSTEM={elementtype.split()}, #_atoms={n_atoms} b={b}", end='\t\n')
	
			
	
	
