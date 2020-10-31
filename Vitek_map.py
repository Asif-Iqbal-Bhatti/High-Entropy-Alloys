#!/usr/bin/env python3
#-----------------------------------------------------------------------------------------
# AUTHOR:: ASIF IQBAL
# USAGE :: python3 sys.argv[0] <perf/final>
# DATED :: 31/10/2020
# NB :: To plot Differential displacement map. ADAPTATION FROM CLOUET and ATOMMAN et al.,
# https://en.wikipedia.org/wiki/Fractional_coordinates#Conversion_to_cartesian_coordinates
# http://www.ruppweb.org/Xray/tutorial/Coordinate%20system%20transformation.htm
#-----------------------------------------------------------------------------------------

import numpy as np
import os, sys, math
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms
from ase import Atoms
from ase.io import write, read
from ase.io.vasp import write_vasp, read_vasp

Struct = sys.argv[1] # Define the reference structure
Per_Z  = True
Per_Y  = False
Per_XY = False
DEBUG  = False

def read_POSCAR_perfect():
	P_pos = []; sum = 0; P_LV1 = []; P_LV2 = []; P_LV3 = []
	reffile = open('POSCAR_perfect','r')
	P_comment  = reffile.readline() # IGNORE first line comment
	P_alat     = float( reffile.readline() )# scale
	P1 = reffile.readline().split()
	P2 = reffile.readline().split()
	P3 = reffile.readline().split() 
	for ai in P1: P_LV1.append(float(ai))
	for bi in P2: P_LV2.append(float(bi))
	for ci in P3: P_LV3.append(float(ci))	
	P_elemtype   = reffile.readline()
	P_atomtypes  = reffile.readline()
	P_Coordtype  = reffile.readline().split()
	P_nat = [int(i) for i in P_atomtypes.split()]
	for i in P_nat: sum = sum + i; P_atoms = sum				
	for x in range(int(P_atoms)):	
		P_pos.append( [ str(i) for i in reffile.readline().rstrip(" ").split()[0:3] ] )
	P_pos = [ [float(P_pos[j][i]) for i in range(3)] for j in range(P_atoms) ]
	reffile.close() 
	return P_atoms,P_pos,P_comment,P_alat,P_LV1,P_LV2,P_LV3,P_elemtype,P_atomtypes,P_Coordtype

# THIS function is redundant
def read_CONTCAR():
	C_pos = []; sum = 0; C_LV1 = []; C_LV2 = []; C_LV3 = []
	inpfile = open('CONTCAR','r')
	C_comment  = inpfile.readline() # IGNORE first line comment
	C_alat       = float( inpfile.readline() )# scale
	C1 = inpfile.readline().split();
	C2 = inpfile.readline().split();
	C3 = inpfile.readline().split();
	for ai in C1: C_LV1.append(float(ai))
	for bi in C2: C_LV2.append(float(bi))
	for ci in C3: C_LV3.append(float(ci))
	C_elemtype   = inpfile.readline()
	C_atomtypes  = inpfile.readline()
	C_Coordtype  = inpfile.readline().split()
	C_nat = [int(i) for i in C_atomtypes.split()]
	for i in C_nat: sum = sum + i; C_atoms = sum				
	for x in range(int(C_atoms)):	
		C_pos.append( [ str(i) for i in inpfile.readline().rstrip(" ").split()[0:3] ] )
	C_pos = [ [float(C_pos[j][i]) for i in range(3)] for j in range(C_atoms) ]		
	inpfile.close() 
	return C_atoms,C_pos,C_comment,C_alat,C_LV1,C_LV2,C_LV3,C_elemtype,C_atomtypes,C_Coordtype
	
#### math.sin function takes argument in radians ONLY
#### Ordering of returning angles variables does matter
def Conversion_2_cartesian_coordinates(LV1,LV2,LV3):
	a=[]; b=[]; c=[];
	for ai in LV1: a.append(float(ai)); x = np.linalg.norm(a)
	for bi in LV2: b.append(float(bi)); y = np.linalg.norm(b)
	for ci in LV3: c.append(float(ci)); z = np.linalg.norm(c)
	P_gamma = math.acos(np.dot(a,b) / (x * y)); # gamma = Cos-1( (a.b)/||a||.||b|| ) 
	P_alpha = math.acos(np.dot(b,c) / (y * z)); # alpha = Cos-1( (b.c)/||b||.||c|| )
	P_beta  = math.acos(np.dot(a,c) / (x * z)); # beta  = Cos-1( (a.c)/||a||.||c|| )
	V = x*y*z * ( np.sqrt(1 + 2 * math.cos(P_alpha) * math.cos(P_beta) * \
	math.cos(P_gamma) - math.cos(P_alpha)**2 - math.cos(P_beta)**2 - math.cos(P_gamma)**2) )	
	M = [ [x, y * math.cos(P_gamma), z * math.cos(P_beta)],
			[0, y * math.sin(P_gamma), z * ( math.cos(P_alpha) - math.cos(P_beta) * math.cos(P_gamma) )/math.sin(P_gamma) ],  
			[0, 0, V / (x * y * math.sin(P_gamma)) ] ]
	return M

#----------------------!!! WRITING TO A DISPLACEMENT FILE !!! ----------------------
def Write_to_file(P_atoms,P_LV1,P_LV2,P_LV3,C_LV1,C_LV2,C_LV3):
	outfile = open("displacement.yaml", 'w')
	outfile.write("{:3d}\n".format(P_atoms) )
	outfile.write("Atom Types refFile Xp Yp Zp displacement Xc-Xp, Yc-Yp, Zc-Zp\n".format() )
	
	Mp = Conversion_2_cartesian_coordinates(P_LV1,P_LV2,P_LV3)
	# This is a technical bcz the lattice vectors of the CONTCAR file does not
	# conform the matrix M that is generated from the "Conversion_..." routine
	# So what I have done here is just directly take the lattice vectors from the 
	# CONTCAR file. I THINK there musr be another way to do that.
	Mc = np.transpose( [C_LV1,C_LV2,C_LV3] )
	
	dict = {}; g = 0;
	for j in range( len( P_elemtype.split() )):
		dict[P_elemtype.split()[j]] =  P_atomtypes.split()[j]; 
	print (dict)
	
	for l in dict:
		for k in range( int(dict[l]) ):
			Xp = P_pos[g][0]; Yp = P_pos[g][1]; Zp = P_pos[g][2]
			Xc = C_pos[g][0]; Yc = C_pos[g][1]; Zc = C_pos[g][2]
			Dp = [Xp, Yp, Zp]; Dc = [Xc, Yc, Zc]
			
			if (P_Coordtype[0] == "Direct" or P_Coordtype[0] == "direct"):
			
				# Conversion to Cartesian coordinates
				Sp = Mp @ np.transpose(Dp); Sc = Mc @ np.transpose(Dc); 
				
				if (Struct == 'perf'):
					#print ("{:2s} {:20.16f} {:20.16f} {:20.16f} {:20.16f} {:20.16f} {:20.16f}".format(l, Sp[0], Sp[1], Sp[2], Sc[0]-Sp[0], Sc[1]-Sp[1], Sc[2]-Sp[2]  )) 			
					outfile.write("{:2s} {:20.16f} {:20.16f} {:20.16f} {:20.16f} {:20.16f} {:20.16f}\n". \
					format(l, Sp[0], Sp[1], Sp[2], Sc[0]-Sp[0], Sc[1]-Sp[1], Sc[2]-Sp[2]  ))
					
				elif (Struct == 'final'):
					#print ("{:2s} {:20.16f} {:20.16f} {:20.16f} {:20.16f} {:20.16f} {:20.16f}".format(l, Sc[0], Sc[1], Sc[2], Sc[0]-Sp[0], Sc[1]-Sp[1], Sc[2]-Sp[2] ))			
					outfile.write("{:2s} {:20.16f} {:20.16f} {:20.16f} {:20.16f} {:20.16f} {:20.16f}\n". \
					format(l, Sc[0], Sc[1], Sc[2], Sc[0]-Sp[0], Sc[1]-Sp[1], Sc[2]-Sp[2] ))	
			
			if  (P_Coordtype[0] == "Cartesian" or P_Coordtype[0] == "cartesian"):
				
				if (Struct == 'perf'):
					#print ("{:15.12f} {:15.12f} {:15.12f} {:15.12f} {:15.12f} {:15.12f}".format(Xp,Yp,Zp, Xc-Xp, Yc-Yp, Zc-Zp ) )			
					outfile.write("{:15.12f} {:15.12f} {:15.12f} {:15.12f} {:15.12f} {:15.12f}\n". \
					format(Xp, Yp, Zp, Xc-Xp, Yc-Yp, Zc-Zp ) )
					
				elif (Struct == 'final'):
					#print ("{:15.12f} {:15.12f} {:15.12f} {:15.12f} {:15.12f} {:15.12f}".format(Xc,Yc,Zc, Xc-Xp, Yc-Yp, Zc-Zp ) )			
					outfile.write("{:15.12f} {:15.12f} {:15.12f} {:15.12f} {:15.12f} {:15.12f}\n". \
					format(Xc, Yc, Zc, Xc-Xp, Yc-Yp, Zc-Zp  ) )	
			g += 1
	outfile.close()
	return 0

def plot_supercell():
	P_ref = read_vasp('POSCAR_perfect')
	P_pos = P_ref.get_positions()
	P_ucell = P_ref.get_cell()
	P_type = P_ref.get_chemical_symbols()
	P_1 = Atoms(positions=P_pos, symbols=P_type, cell=P_ucell, pbc=[True,True,True])
	
	C_ref = read_vasp('CONTCAR')
	C_pos = C_ref.get_positions()
	C_ucell = C_ref.get_cell()
	C_type = C_ref.get_chemical_symbols()
	C_1 = Atoms(positions=C_pos, symbols=C_type, cell=C_ucell, pbc=[True,True,True])

	fig, ax = plt.subplots(2, 1, figsize=(5, 5) )
	plot_atoms(P_1, ax[0], radii=0.3, rotation=('0x,0y,00z'))
	ax[0].set_title("POSCAR perfect")
	plot_atoms(C_1, ax[1], radii=0.3, rotation=('0x,0y,00z'))
	ax[1].set_title("CONTCAR")	
	plt.show()

#----------------------!!! FOR BCC STRUCTURE ONLY
def Vitek_map(P_atoms): 
	alat = 3.38 # Angstrom
	bz = alat/2 * np.sqrt(3)
	# Half distance between 1st and 2nd n.n.
	Rcut2 = ( 1/2 * alat * (np.sqrt(6)/3 + np.sqrt(2)) )**2
	print ("{} {}".format(alat, bz) )
	P_x = []; P_y = []; P_z = []; tmp = []
	P_gradx = []; P_grady = []; P_gradz = []; 
	
	disp_file = open('displacement.yaml', 'r')
	disp_file.readline().strip("\n").split() # first line comment
	disp_file.readline().strip("\n").split() # 2nd line comment
	for j in range(P_atoms):
		tmp.append( disp_file.readline().strip("\n").split() )
	disp_file.close()

	for i in range(0,P_atoms,1):
		P_x.append(tmp[i][1] ); 
		P_y.append(tmp[i][2] )
		P_z.append(tmp[i][3] )
		P_gradx.append(tmp[i][4] )
		P_grady.append(tmp[i][5] )
		P_gradz.append(tmp[i][6] )
	P_x = np.array(P_x, dtype=np.float64)	
	P_y = np.array(P_y, dtype=np.float64)	
	P_z = np.array(P_z, dtype=np.float64)	
	P_gradx = np.array(P_gradx, dtype=np.float64)	
	P_grady = np.array(P_grady, dtype=np.float64)	
	P_gradz = np.array(P_gradz, dtype=np.float64)	
	
	# Dimension of the simulation box
	x0 = min(P_x)
	x1 = max(P_x)	
	y0 = min(P_y)
	y1 = max(P_y)
	x_diff = x1-x0
	y_diff = y1-y0	
	
	### Periodic boundary condition along Z
	if (Per_Z):
		for i in range (P_atoms):
			P_z[i] = P_z[i]  - math.floor( P_z[i] /bz )*bz

  ### File containing points of the different planes 
	points = open("points.yaml","w")
	z0 = min(P_z) + bz/6.
	
	### Atom positions corresponding to 1st (111) plane
	### For any layers number of layers in the Z direction the following definition
	### could be changed. FOr example in case of HEAs where the Z direction is twice
	### the <111>.
	for i in range (P_atoms):
		if ( 3.0*np.mod(P_z[i]-z0, bz) < bz ):
			points.write('{:14.6f} {:14.6f} {:14.6f}\n'.format( P_x[i], P_y[i], P_z[i] ) )
	points.write("\n")

	### Atom positions corresponding to 2nd (111) plane
	for i in range (P_atoms):
		if ( (3.0*np.mod(P_z[i]-z0,bz) >= bz) and (3.0*np.mod(P_z[i]-z0,bz) < 2.0*bz) ):
			points.write('{:14.6f} {:14.6f} {:14.6f}\n'.format( P_x[i], P_y[i], P_z[i] ) )
	points.write("\n")
	
	### Atom positions corresponding to 3rd (111) plane
	for i in range (P_atoms):
		if (3.0*np.mod(P_z[i]-z0,bz) >= 2.0*bz):
			points.write('{:14.6f} {:14.6f} {:14.6f}\n'.format( P_x[i], P_y[i], P_z[i] ) )
	
	points.close()
	
	### Output file for screw components 
	res = open("res.yaml","w")
	res.write('# [0.0, 0.0, {}] (Burgers vector)\n'.format(bz)         )
	res.write('# [0.0, 0.0, 1.0] (line direction)\n'           )
	res.write('# 1-3: RefFile x, y, and z coordinates\n'             )
	res.write('# 4-6: FinalFile x, y, and z coordinates\n'           )
	res.write('#   7: displacement difference uz2-uz1 along b\n' )
	res.write('# 8-9: atom indexes i and j\n'                    )
	
	### Computing the arrows in the <111> directions
	z0 = min( P_z ) + 2.0*bz
	dRij_shift  = 0.0   # Shift to apply to position of j atom for PBC
	imVitek = P_atoms   # Number of atoms for Vitek maps including periodic images
	at = [P_LV1, P_LV2, P_LV3]
	inv_at = np.linalg.inv( at ) # INVERSE of a lattice vector

	for i in range( P_atoms-1 ):  # i>j
		if ( P_z[i] > z0): continue
		for j in range(i+1, P_atoms):
			if ( P_z[j] > z0): continue
			dRij = [ P_x[j] - P_x[i], P_y[j] - P_y[i], P_z[j] - P_z[i] ]

			if (Per_XY):
				dSij = np.dot( inv_at, dRij )
				if ( ( round( dSij[0] ) != 0 ) or ( round(dSij[1]) != 0 ) ):
					if (DEBUG):
						res.write('before: dRij = {}'.format( dRij ) )
						dSij = dSij - round( dSij )
						res.write('after: dRij = {}'.format( np.dot( at, dSij )) )
						res.write('=> dR2 = {}', dRij[0]**2 + dRij[1]**2 )
						res.write("")
					dSij = dSij - round( dSij )
					dRij_shift = np.dot( at, dSij ) - dRij
					dRij = dRij + dRij_shift
					Per_Y = True # We are considering a periodic image of j
				else:
					dRij_shift = 0.0
					Per_Y = False # We are not considering a periodic image of j

			R2 = dRij[0]**2 + dRij[1]**2
			if (R2 > Rcut2): continue
			if (R2 <= 1.0E-8): continue
			if (Per_XY):
				imVitek = imVitek + 1
				jVitek = imVitek
			else:
				jVitek = j
			inv_R = 1/np.sqrt(R2)
			Dgradr = [ P_gradx[j] - P_gradx[i], P_grady[j] - P_grady[i] ]  
			Dgradz = P_gradz[j] - P_gradz[i]
			Dgradz = Dgradz - bz*round(Dgradz/bz)
			
			res.write("{:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:2d} {:2d}\n". \
			format(P_x[i], P_y[i], P_z[i], P_x[j], P_y[j], P_z[j] + dRij_shift, Dgradz, i, jVitek ) )
			
			res.write("{:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:12.6f} {:2d} {:2d}\n". \
			format(P_x[i], P_y[i], P_z[i], P_x[j], P_y[j], P_z[j] + dRij_shift, Dgradr, i, jVitek ) )			
	res.close()
	
#---------------------------------- MAIN ENGINE ------------------------------------

if __name__ == "__main__":
	P_atoms,P_pos,P_comment,P_alat,P_LV1,P_LV2,P_LV3,P_elemtype,P_atomtypes,P_Coordtype = read_POSCAR_perfect();	
	C_atoms,C_pos,C_comment,C_alat,C_LV1,C_LV2,C_LV3,C_elemtype,C_atomtypes,C_Coordtype = read_CONTCAR();	

	### CHECKING IF THE ATOMS AND ELEMENTS ARE EQUAL in both files
	if (P_atoms == C_atoms):
		print("atoms are equal")
	else:
		exit("atoms are not equal")	
	if (P_elemtype.split()  == C_elemtype.split() ):
		print ("Elements are equal")
	else:
		exit("Elements are not equal")
	if ( P_Coordtype[0] == C_Coordtype[0] ):
		print ("Coordinates types are equal")
	else:
		exit("Coordinates types are not equal")

	Write_to_file(P_atoms,P_LV1,P_LV2,P_LV3,C_LV1,C_LV2,C_LV3)
	Vitek_map(P_atoms)
	
	
	
	
	
	
	
