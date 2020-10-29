#!/usr/bin/env python3

#---------------------------------------------------------------------------------
''' 
# AUTHOR::: Asif Iqbal
# DATED ::: 28/10/2020
# GITHUB::: @asif_em
# USAGE		:: python sys.argv[0] bestsqs.out <c/d> 
# This script reads bestsqs.out generated from ATAT mcsqs code.
# Please visit https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/
'''
#---------------------------------------------------------------------------------

import numpy as np
import os, sys
from numpy import linalg as LA
from termcolor import colored

def load_bestsqs():
	coord=[]; pos=[]; elem=[]; elem_O=[]; atm_num=[]; count=0; CX=[]; LV=[]
	
	content = len( open(sys.argv[1]).readlines( ) )
	atoms = content-6 # substract 6 lines because first six lines are vectors
	print("{:_^60s}".format("*"))
	print("# of atoms {} in the system".format(atoms) )

	fdata1 = open(sys.argv[1],'r')
	#---
	Cart_Vec1 = fdata1.readline().split(); Cart_VX = [float(i) for i in Cart_Vec1]
	Cart_Vec2 = fdata1.readline().split(); Cart_VY = [float(i) for i in Cart_Vec2]
	Cart_Vec3 = fdata1.readline().split(); Cart_VZ = [float(i) for i in Cart_Vec3]
	Lat_Vec1 	= fdata1.readline().split(); Lat_VX = [float(i) for i in Lat_Vec1]
	Lat_Vec2 	= fdata1.readline().split(); Lat_VY = [float(i) for i in Lat_Vec2]
	Lat_Vec3 	= fdata1.readline().split(); Lat_VZ = [float(i) for i in Lat_Vec3]	
	CX.append(Cart_VX); CX.append(Cart_VY); CX.append(Cart_VZ)
	LV.append(Lat_VX); LV.append(Lat_VY); LV.append(Lat_VZ)		
	print ("Cartesian Vector: {}".format(CX))
	print ("Lattice Vector: {}".format(LV))
	POS_vec = np.dot(CX, LV)
	
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
	fdata1.close()
	
	elem_atom = zip(elem, atm_num)
	dict1 = {a : b for a,b in elem_atom}
	print ("{}".format( dict1 ) )
	
#------------------- Converting bestsqs.out to POSCAR file (VASP5)
	fdata2 = open('POSCAR_bestsqs','w')
	
	fdata2.write("POSCAR generated from bestsqs\n")
	fdata2.write("{:6.6f}\n".format(1.0))
	for j in range(3):
		fdata2.write("{:12.9f} {:12.9f} {:12.9f}\n".format(POS_vec[j][0], POS_vec[j][1], POS_vec[j][2]))
	for i in elem_O: fdata2.write("{:6s}".format(i)); 
	fdata2.write("\n")	
	for i in atm_num: fdata2.write("{:4d}".format(i)); 
	fdata2.write("\n")	

#------------------------- In Cartesian UNITS ----------------#

	if (sys.argv[2] == 'c' or sys.argv[2] == 'C' or sys.argv[2] == 'Cartesian'):
		fdata2.write("{:12.9s}\n".format("Cartesian"))

		for i in range(atoms):
			#a1 = pos[i][0]*CX[0][0] + pos[i][1]*CX[1][0] + pos[i][2]*CX[2][0]
			#a2 = pos[i][0]*CX[0][1] + pos[i][1]*CX[1][1] + pos[i][2]*CX[2][1]
			#a3 = pos[i][0]*CX[0][2] + pos[i][1]*CX[1][2] + pos[i][2]*CX[2][2]
			A = np.dot( pos[i][:], CX )
			fdata2.write("{:15.12f} {:15.12f} {:15.12f}\n".format( A[0], A[1], A[2] ) )
		
#------------------------- In fractional/reduced UNITS ----------------#	

	elif (sys.argv[2] == 'd' or sys.argv[2] == 'D' or sys.argv[2] == 'Direct'):	
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

	fdata2.close()
	print (colored("File has been generated in {}".format(sys.argv[2]), "yellow" ))
	
if __name__ == "__main__":
	load_bestsqs()
	
