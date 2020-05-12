#!/usr/bin/env python3
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# USAGE :: ./python3 <alat> <POSCAR/CONTCAR>
# <POSCAR/CONTCAR> should be in Cartesian coordinates
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

import numpy as np
import os, sys

scale = sys.argv[1]
alat = float(scale)
POSCAR = sys.argv[2]
pos = []; kk = []; lattice = []; sum = 0
print('<POSCAR/CONTCAR> should be in Cartesian coordinates: \n')
print('Reading File ... : \n')

file = open(POSCAR,'r')

firstline   = file.readline() # IGNORE first line comment
secondline = file.readline() # scale
Latvec1 = file.readline().split()
print(float(Latvec1[0])/alat,float(Latvec1[1])/alat,float(Latvec1[2])/alat)
Latvec2 = file.readline().split()
print(float(Latvec2[0])/alat,float(Latvec2[1])/alat,float(Latvec2[2])/alat)
Latvec3 = file.readline().split()
print(float(Latvec3[0])/alat,float(Latvec3[1])/alat,float(Latvec3[2])/alat)
elementtype=file.readline()
print ("Types of elements:", str(elementtype), end = '\n')
numofatoms=file.readline()
Coordtype=file.readline()
print ("Coordtype:", (Coordtype), end = '\n')	

print (">>>>>>>>>------------------- add # of Atoms--------------------")
nat = numofatoms.split()
nat = [int(i) for i in nat]
print (nat)
for i in nat:
	sum = sum + i
numberofatoms = sum
print ("Number of atoms:", (numberofatoms), end = '\n')

#>>>>>>>>>---------------Atomic positions------------------")				
for x in range(int(numberofatoms)):
	coord = file.readline().split()
	coord = [float(i) for i in coord]
	pos = pos + [coord]
pos = np.array(pos)
file.close()	
#>>>>>>>>>---------------Atomic positions------------------")

fdata2 = open('POSCAR_scale','w')
fdata2.write(firstline)
fdata2.write("{:12.6f}\n".format(alat) )
fdata2.write("{:15.12f} {:15.12f} {:15.12f}\n".format( float(Latvec1[0])/alat,float(Latvec1[1])/alat,float(Latvec1[2])/alat) )
fdata2.write("{:15.12f} {:15.12f} {:15.12f}\n".format( float(Latvec2[0])/alat,float(Latvec2[1])/alat,float(Latvec2[2])/alat) )
fdata2.write("{:15.12f} {:15.12f} {:15.12f}\n".format( float(Latvec3[0])/alat,float(Latvec3[1])/alat,float(Latvec3[2])/alat) )
fdata2.write(elementtype)
fdata2.write(numofatoms)
fdata2.write(Coordtype)
for x in range(int(numberofatoms)):
	fdata2.write("{:12.9f} {:12.9f} {:12.9f}\n".format(pos[x][0]/alat, pos[x][1]/alat, pos[x][2]/alat) )	
fdata2.close()	
#---------------------------------------------------------
