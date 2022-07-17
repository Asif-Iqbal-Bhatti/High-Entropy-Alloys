#!/usr/bin/env python3

#---------------------------------------------------------------------------------
''' 
# AUTHOR :: Asif Iqbal
# DATED  :: 28/10/2020
# GITHUB :: @aib_em
# USAGE	 :: python sys.argv[0] bestsqs.out <c/d, Cartesian or Direct> 
# This script reads bestsqs.out generated from ATAT mcsqs code.
# Visit https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/
'''
#---------------------------------------------------------------------------------

import numpy as np
import os, sys
from numpy import linalg as LA

def load_bestsqs():
	coord=[]
	pos=[]
	elem=[]
	elem_O=[]
	atm_num=[]
	count=0
	CX=[]
	LV=[]

	content = len( open(sys.argv[1]).readlines( ) )
	atoms = content-6 # substract 6 lines because first six lines are vectors
	print("{:_^60s}".format("*"))
	print(f"# of atoms {atoms} in the system")

	with open(sys.argv[1],'r') as fdata1:
		CX.extend(
			(
				[float(i) for i in fdata1.readline().split()],
				[float(i) for i in fdata1.readline().split()],
				[float(i) for i in fdata1.readline().split()],
			)
		)

		LV.extend(
			(
				[float(i) for i in fdata1.readline().split()],
				[float(i) for i in fdata1.readline().split()],
				[float(i) for i in fdata1.readline().split()],
			)
		)

		POS_vec = np.dot(CX, LV)
		print(POS_vec)

		#--- Sorting coordinates and element types
		for i in range(atoms):
			lines = fdata1.readline().split()
			coord.append(lines)
			elem.append(coord[i][3])
		elem = sorted( list(set(elem)) )

		#--- Counting number of elements and sorting coordinates
		for i in range(len(elem)):
			for j in range(atoms):
				if (coord[j][3] == elem[i]):
					count = count + 1
					pos.append(coord[j])
			elem_O.append(elem[i])
			atm_num.append(count)	
			count = 0 # RESETING COUNTER FOR CORRECT ATOMS NUMBERING
		pos = [ [float(pos[j][i]) for i in range(3)] for j in range(atoms) ]
	elem_atom = zip(elem, atm_num)
	dict1 = dict(elem_atom)
	print(f"{dict1}")

	with open('POSCAR_bestsqs','w') as fdata2:
		fdata2.write("POSCAR generated from bestsqs\n")
		fdata2.write("{:6.6f}\n".format(1.0))
		for j in range(3):
			fdata2.write("{:12.9f} {:12.9f} {:12.9f}\n".format(POS_vec[j][0], POS_vec[j][1], POS_vec[j][2]))
		for i in elem_O: fdata2.write("{:6s}".format(i));
		fdata2.write("\n")
		for i in atm_num: fdata2.write("{:4d}".format(i));
		fdata2.write("\n")	

			#------------------------- In Cartesian UNITS ----------------#

		if sys.argv[2] in ['c', 'C', 'Cartesian']:
			fdata2.write("{:12.9s}\n".format("Cartesian"))

			for i in range(atoms):
				#a1 = pos[i][0]*CX[0][0] + pos[i][1]*CX[1][0] + pos[i][2]*CX[2][0]
				#a2 = pos[i][0]*CX[0][1] + pos[i][1]*CX[1][1] + pos[i][2]*CX[2][1]
				#a3 = pos[i][0]*CX[0][2] + pos[i][1]*CX[1][2] + pos[i][2]*CX[2][2]
				A = np.dot( pos[i][:], CX )
				fdata2.write("{:15.12f} {:15.12f} {:15.12f}\n".format( A[0], A[1], A[2] ) )

		elif sys.argv[2] in ['d', 'D', 'Direct']:	
			fdata2.write("{:12.9s}\n".format("Direct"))
			u = np.cross(POS_vec[1], POS_vec[2])
			v = np.cross(POS_vec[0], POS_vec[2])
			w = np.cross(POS_vec[0], POS_vec[1])
			V = np.array([ POS_vec[0],POS_vec[1],POS_vec[2] ] )
			print ("Supercell Volume :: {:9.6f} & Volume/atom :: {:9.6f}".format(LA.det(V), LA.det(V)/atoms ) )
			Vx = float( np.dot(POS_vec[0],u) )
			Vy = float( np.dot(POS_vec[1],v) ) 
			Vz = float( np.dot(POS_vec[2],w) )			

			for i in range(atoms):
				#a1 = pos[i][0]*CX[0][0] + pos[i][1]*CX[1][0] + pos[i][2]*CX[2][0]
				#a2 = pos[i][0]*CX[0][1] + pos[i][1]*CX[1][1] + pos[i][2]*CX[2][1]
				#a3 = pos[i][0]*CX[0][2] + pos[i][1]*CX[1][2] + pos[i][2]*CX[2][2]
				A = np.dot( pos[i][:], CX )
				fdata2.write("{:12.9f} {:12.9f} {:12.9f}\n".format(np.dot(A,u)/Vx, np.dot(A,v)/Vy, np.dot(A,w)/Vz ) )	

	print(f"File has been generated in {sys.argv[2]} coord")
	
if __name__ == "__main__":
	load_bestsqs()
	
