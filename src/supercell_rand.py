#!/usr/bin/env python3

############################################################################
# USAGE  :: python3 sys.argv[0] bcc/fcc c/d
# Author :: Asif Iqbal
# DATED  :: 29/10/2020
############################################################################

import numpy as np
import sys, random
from numpy import linalg as LA

#------------------------------INPUT PARAMETRS------------------------------------#
lattice_parameter = 3.405 
supercellx = 5
supercelly = 5
supercellz = 5
#----------------------------------------------------------------------------------#
def HEAs_supercell():
	cartesian_units=[]
	count=0
	lattice_vector = np.array([[1,0,0],
														[0,1,0],
														[0,0,1]])*lattice_parameter

	lattice_bcc = np.array([[0,0,0],
												[0.5,0.5,0.5]])*lattice_parameter

	lattice_fcc = np.array([[0,0,0],
												[0.0,0.5,0.5],											
												[0.5,0.0,0.5],											
												[0.5,0.5,0.0]])*lattice_parameter
	fcc = lattice_fcc;
	if (sys.argv[1]=='fcc'): b = lattice_parameter*np.sqrt(2)/2.0; print("Burgers vector::", b)
	bcc = lattice_bcc
	if (sys.argv[1]=='bcc'): b = lattice_parameter*np.sqrt(3)/2.0; print("Burgers vector::", b)

	for i in range(supercellx):										
		for j in range(supercelly):
			for k in range(supercellz):
				atom_position = np.array([i,j,k])
				cartesian_basis = np.inner(lattice_vector.T, atom_position)
				if (sys.argv[1] == 'bcc'):
					for atom in lattice_bcc:
						cartesian_units.append(cartesian_basis + atom)
						count+=1
				elif (sys.argv[1] == 'fcc'):		
					for atom in lattice_fcc:
						cartesian_units.append(cartesian_basis + atom)
						count+=1			

	with open("POSCAR","w") as POSCAR:
		POSCAR.write(f'#{sys.argv[1]}\n')
		POSCAR.write('{:6.6f}\n'.format(1.0))
		POSCAR.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lattice_vector[0][0]*supercellx,lattice_vector[0][1]*supercellx,lattice_vector[0][2]*supercellx ))
		POSCAR.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lattice_vector[1][0]*supercelly,lattice_vector[1][1]*supercelly,lattice_vector[1][2]*supercelly ))
		POSCAR.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lattice_vector[2][0]*supercellz,lattice_vector[2][1]*supercellz,lattice_vector[2][2]*supercellz ))
		print("{:12.9f} {:12.9f} {:12.9f}".format(lattice_vector[0][0]*supercellx,lattice_vector[0][1]*supercellx,lattice_vector[0][2]*supercellx ))
		print("{:12.9f} {:12.9f} {:12.9f}".format(lattice_vector[1][0]*supercelly,lattice_vector[1][1]*supercelly,lattice_vector[1][2]*supercelly ))
		print("{:12.9f} {:12.9f} {:12.9f}".format(lattice_vector[2][0]*supercellz,lattice_vector[2][1]*supercellz,lattice_vector[2][2]*supercellz ))
		POSCAR.write('Ta\n')
		POSCAR.write(f'{count}\n')
	#------------------------- In Cartesian UNITS -------------------------#
		if sys.argv[2] in ["c", "C"]:
			POSCAR.write("Cartesian\n")
			for cartesian_unit in cartesian_units:
			#print("{:12.9f} {:12.9f} {:12.9f}".format(cartesian_units[i][0], cartesian_units[i][1], cartesian_units[i][2]))
				POSCAR.write(
					"{:12.9f} {:12.9f} {:12.9f}\n".format(
						cartesian_unit[0], cartesian_unit[1], cartesian_unit[2]
					)
				)
	#------------------------- In fractional/Direct UNITS -------------------------#		
		u = np.cross(lattice_vector[1]*supercelly, lattice_vector[2]*supercellz)
		v = np.cross(lattice_vector[0]*supercellx, lattice_vector[2]*supercellz)
		w = np.cross(lattice_vector[0]*supercellx, lattice_vector[1]*supercelly)
		V = np.array([ lattice_vector[0]*supercelly,lattice_vector[1]*supercelly,lattice_vector[2]*supercellz ] )
		print ("Volume of the cell::", LA.det(V) )
		Vx = np.inner(lattice_vector[0]*supercellx,u)
		Vy = np.inner(lattice_vector[1]*supercelly,v)
		Vz = np.inner(lattice_vector[2]*supercellz,w)

		if sys.argv[2] in ["d", "D"]:
			POSCAR.write("Direct\n")
			for cartesian_unit_ in cartesian_units:
				POSCAR.write(
					"{:12.9f} {:12.9f} {:12.9f}\n".format(
						np.dot(cartesian_unit_, u) / Vx,
						np.dot(cartesian_unit_, v) / Vy,
						np.dot(cartesian_unit_, w) / Vz,
					)
				)


	#------------------------- Randomly distribute atoms for HEA -------------------------#
	with open("newPOSCAR","w") as fdata:
		for _ in range(int(1e3)):
			randomArrx = random.sample(range(count), count)
		if (count%5 == 0): 	
			fdata.write(f'#{sys.argv[1]}\n')
			fdata.write('{:6.6f}\n'.format(1.0))
			fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lattice_vector[0][0]*supercellx,lattice_vector[0][1]*supercellx,lattice_vector[0][2]*supercellx ))
			fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lattice_vector[1][0]*supercelly,lattice_vector[1][1]*supercelly,lattice_vector[1][2]*supercelly ))
			fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lattice_vector[2][0]*supercellz,lattice_vector[2][1]*supercellz,lattice_vector[2][2]*supercellz ))
			fdata.write('A B C D E\n')
			fdata.write('{0} {0} {0} {0} {0}\n'.format((int(count/5))) )

			#------------------------- In Cartesian UNITS -------------------------#
			if sys.argv[2] in ["c", "C"]:
				fdata.write("Cartesian\n")
				print(
					f'# of atoms per element -> {count / 5}. File generated in Cartesian coordinates'
				)

				for i in range(len(cartesian_units) ):
					#print("{:12.9f} {:12.9f} {:12.9f}".format(cartesian_units[randomArrx[i]][0],cartesian_units[randomArry[i]][1],cartesian_units[randomArrz[i]][2] ) )
					fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(cartesian_units[randomArrx[i]][0],cartesian_units[randomArrx[i]][1],cartesian_units[randomArrx[i]][2] ) )

			#------------------------- In fractional/reduced UNITS ----------------#
			if sys.argv[2] in ["d", "D"]:
				fdata.write("Direct\n")
				print(
					f'# of atoms per element -> {count / 5}. File generated in Direct coordinates'
				)

				for i in range(len(cartesian_units) ):
					fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(np.dot(cartesian_units[randomArrx[i]],u)/Vx, np.dot(cartesian_units[randomArrx[i]],v)/Vy, np.dot(cartesian_units[randomArrx[i]],w)/Vz ) )
		else: 
			print(
				f'{count / 5} is not an integer number for equal composition -> HEAs not generated'
			)		
	
def help():
	print('A simple script to generate FCC or BCC supercell for HEAs.')
	print('To execute just run python3 sys.argv[0] <bcc/fcc> <c/d>.')
	print('THIS script is valid for equimolar composition !!!')	
	print('HEAs consists of five or more elements. The elements has already been typed into the')
	print('script just change according to your need also lattice vectors should be')
	print('equal and atomic composition should corresponds to integer multiple of atoms.')
	
if __name__ == '__main__':
	if len(sys.argv) < 3 or len(sys.argv) > 3:
		help()
	else:
		HEAs_supercell()
			
			
			
			
			
			
			
			
