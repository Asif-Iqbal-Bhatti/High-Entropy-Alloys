#!/usr/bin/env python3

''' 
>>> This script reads bestsqs.out generated from ATAT mcsqs code.
>>> Please visit https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/
>>> Authors :: Asif Iqbal 
>>> DATED 	:: 09/05/2020 PARIS, FRANCE
>>> USAGE		:: python sys.argv[0] bestsqs.out <c/d> 
>>> 
 '''
 
import numpy as np
import os, sys
from numpy import linalg as LA

def load_bestsqs():
	coord=[]; pos=[]; elem=[]; elem_O=[]; atm_num=[]; count=0; CX=[]; LV=[]
	
	content = len( open(sys.argv[1]).readlines( ) )
	atoms = content-6 # substract 6 because first sixe lines are vectors
	print ("-"*50)
	print("# of atoms", atoms)

	fdata1 = open(sys.argv[1],'r')
	#---
	Cart_Vec1 = fdata1.readline().split(); Cart_VX = [float(i) for i in Cart_Vec1]; print (Cart_VX)
	Cart_Vec2 = fdata1.readline().split(); Cart_VY = [float(i) for i in Cart_Vec2]; print (Cart_VY)
	Cart_Vec3 = fdata1.readline().split(); Cart_VZ = [float(i) for i in Cart_Vec3]; print (Cart_VZ)
	CX.append(Cart_VX); CX.append(Cart_VY); CX.append(Cart_VZ)
	#---
	Lat_Vec1 = fdata1.readline().split(); Lat_VX = [float(i) for i in Lat_Vec1]; print (Lat_VX)
	Lat_Vec2 = fdata1.readline().split(); Lat_VY = [float(i) for i in Lat_Vec2]; print (Lat_VY)
	Lat_Vec3 = fdata1.readline().split(); Lat_VZ = [float(i) for i in Lat_Vec3]; print (Lat_VZ)
	LV.append(Lat_VX); LV.append(Lat_VY); LV.append(Lat_VZ)	
	
	#--- Sorting coordinates and element types
	for i in range(atoms):
		lines = fdata1.readline().split()
		coord.append(lines)
		elem.append(coord[i][3])
	elem = list(set(elem))
	print (elem)
	
	#--- Counting number of elements and sorting
	for i in range(len(elem)):
		for j in range(atoms):
			if (coord[j][3] == elem[i]):
				count = count + 1
				pos.append(coord[j])
		elem_O.append(elem[i])
		atm_num.append(count)	
		count = 0 # RESETING COUNTER FOR CORRECT ATOMS NUMBERING
	print (atm_num)
	pos = [ [float(pos[j][i]) for i in range(3)] for j in range(atoms) ]
	fdata1.close()

#--- Converting bestsqs.out to POS=POSCAR output file for VASP5
	fdata2 = open('POSCAR_bestsqs','w')
	
	fdata2.write("POSCAR generated from bestsqs\n")
	fdata2.write("{:6.6f}\n".format(1.0))
	POS_vec = np.dot(CX, LV)
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
			a1 = pos[i][0]*CX[0][0] + pos[i][1]*CX[1][0] + pos[i][2]*CX[2][0]
			a2 = pos[i][0]*CX[0][1] + pos[i][1]*CX[1][1] + pos[i][2]*CX[2][1]
			a3 = pos[i][0]*CX[0][2] + pos[i][1]*CX[1][2] + pos[i][2]*CX[2][2]
			fdata2.write("{:15.12f} {:15.12f} {:15.12f}\n".format( a1, a2, a3 ) )
	else:	
#------------------------- In fractional/reduced UNITS ----------------#
		fdata2.write("{:12.9s}\n".format("Direct"))
		u = np.cross(POS_vec[1], POS_vec[2])
		v = np.cross(POS_vec[0], POS_vec[2])
		w = np.cross(POS_vec[0], POS_vec[1])
		V = np.array([ POS_vec[0],POS_vec[1],POS_vec[2] ] ); print ("Volume of the cell::", LA.det(V) )
		Vx = np.inner(POS_vec[0],u); Vy = np.inner(POS_vec[1],v); Vz = np.inner(POS_vec[2],w)			
		
		for i in range(atoms):
			ax = [float(pos[i][0])*Cart_VX[0], float(pos[i][1])*Cart_VY[0], float(pos[i][2])*Cart_VZ[0] ]
			ay = [float(pos[i][0])*Cart_VX[1], float(pos[i][1])*Cart_VY[1], float(pos[i][2])*Cart_VZ[1] ]
			az = [float(pos[i][0])*Cart_VX[2], float(pos[i][1])*Cart_VY[2], float(pos[i][2])*Cart_VZ[2] ]
			fdata2.write("{:12.9f} {:12.9f} {:12.9f}\n".format(np.dot(ax,u)/Vx, np.dot(ay,v)/Vy, np.dot(az,w)/Vz ) )	

	fdata2.close()
	
if __name__ == "__main__":
	load_bestsqs()
	
