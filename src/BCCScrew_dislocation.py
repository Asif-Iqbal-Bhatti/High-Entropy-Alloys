#!/usr/bin/env python3

'''!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!!! USAGE ::: python3 <fcc/bcc> <c/d> <sup/screw>
!!! AUTHOR::: Asif Iqbal
!!! DATED ::: 12/10/2020
!!! GITHUB::: @asif_em
!!! USE AT YOUR OWN RISK. NOT EVEN IMPLIED WARRANTY WHATSOEVER
!!! CAREFULLY CHECK THE GEOMETRY BEFORE SUBMITTING TO DFT CALCULATION.
!!! This script creates a screw dislocations in a BCC supercell
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------'''

import numpy as np
import sys, random, math
from numpy import linalg as LA

#------------------------------INPUT PARAMETRS------------------------------------#
# Fazakas, E. et al. Experimental and theoretical study of Ti20Zr20Hf 20Nb20X20 
# (X = v or Cr) refractory high-entropy alloys. Int. J. Refract. Met. Hard Mater. 47, 
# 131–138 (2014).
pivalue = 2.0*np.arcsin(1.0) 
lat_par = 3.405 # lattice parameter
Nx = 7; 
Ny = 11;
Nz = 1
#----------------------------------------------------------------------------------- 

def HEAs_supercell():
	SX = 5
	SY = 5
	SZ = 5
	cart_atoms=[]; count=0
	lat_vec = np.array([[1,0,0],
											[0,1,0],
											[0,0,1]])*lat_par

	lattice_bcc = np.array([[0,0,0],
												[0.5,0.5,0.5]])*lat_par

	lattice_fcc = np.array([[0,0,0],
												[0.0,0.5,0.5],											
												[0.5,0.0,0.5],											
												[0.5,0.5,0.0]])*lat_par										
	fcc = lattice_fcc; 
	if (sys.argv[1]=='fcc'): b = lat_par*np.sqrt(2)/2.0; print("burger vector::", b)
	bcc = lattice_bcc; 
	if (sys.argv[1]=='bcc'): b = lat_par*np.sqrt(3)/2.0; print("burger vector::", b)

	for i in range(SX):										
		for j in range(SY):
			for k in range(SZ):
				atom_position = np.array([i,j,k])
				cartesian_basis = np.inner(lat_vec.T, atom_position)
				if (sys.argv[1] == 'bcc'):
					for atom in lattice_bcc:
						cart_atoms.append(cartesian_basis + atom)
						count+=1
				elif (sys.argv[1] == 'fcc'):		
					for atom in lattice_fcc:
						cart_atoms.append(cartesian_basis + atom)
						count+=1
	return cart_atoms, count, lat_vec, SX, SY, SZ		
								
##################################################################################
#------------------ Create Z=<111> surface for Screw Dislocation in bcc ---------#
##################################################################################				

def BCCgen_111_perfect():
	#a1 = [-1-12]; a2=[1-10]; a3=1/2[111]
	Rx = np.sqrt(6.0)        # a*[-1-12]=X
	Ry = np.sqrt(2.0)        # a*[1-10] =Y
	Rz = np.sqrt(3.0)/2.0    # a*[111]  =Z --> a*<111>/2	

	m = 15
	n = 9
	cart_atoms = []
	count=0
	pos_screw = []
	pos_perfect = [];
#	basis = np.array( [ [ -n/3 -1/(6*m), -n/3 -1/(6*m), (2*n)/3 - 1/(6*m) ],
#			   [ -n/6 + m/2 + 1/4 -1/(12*m), -n/6 - m/2 + 1/4 -1/(12*m), (2*n)/6 + 1/4 -1/(12*m) ],
#			   [ 0.5, 0.5, 0.5] ] )*lat_par

	basis = np.array( [ [ Rx, 0, 0 ],
											[ 0, Ry, 0 ],
											[ 0, 0, Rz] ] )*lat_par
# direct coordinates -- > Cartesian	
	basis_atoms = np.array([[0       * Rx ,0   * Ry ,0       *Rz ],
													[1.0/2.0 * Rx ,0.5 * Ry ,0       *Rz ],
													[1.0/3.0 * Rx ,0.0 * Ry ,2.0/3.0 *Rz ],
													[5.0/6.0 * Rx ,0.5 * Ry ,2.0/3.0 *Rz ],
													[1.0/6.0 * Rx ,0.5 * Ry ,1.0/3.0 *Rz ],
													[2.0/3.0 * Rx ,0.0 * Ry ,1.0/3.0 *Rz ]] )*lat_par 

	print ("Original basis vectors :: \n", basis, end="\n\n")

	for i in range(Nx):										
		for j in range(Ny):
			for k in range(Nz):
				atom_position = np.array([i,j,k])
				cartesian_basis = np.inner(basis.T, atom_position)
				for atom in basis_atoms:
					cart_atoms.append(cartesian_basis + atom)
					count+=1

	pos_perfect = cart_atoms
	lx = Nx*Rx*lat_par;
	ly = Ny*Ry*lat_par;
	lz = Nz*Rz*lat_par;

	#------------------------------ Writing to a POSCAR_perfect file ------------------------------
	with open("POSCAR_perfect","w") as fdata2:
		fdata2.write(f'#{sys.argv[1]}\n')
		fdata2.write('{:6.6f}\n'.format(lat_par))
		fdata2.write("{:12.9f} {:12.9f} {:12.9f}\n".format(basis[0][0]*Nx/lat_par,basis[0][1]*Nx/lat_par,basis[0][2]*Nx/lat_par ))
		fdata2.write("{:12.9f} {:12.9f} {:12.9f}\n".format(basis[1][0]*Ny/lat_par,basis[1][1]*Ny/lat_par,basis[1][2]*Ny/lat_par ))
		fdata2.write("{:12.9f} {:12.9f} {:12.9f}\n".format(basis[2][0]*Nz/lat_par,basis[2][1]*Nz/lat_par,basis[2][2]*Nz/lat_par ))
		fdata2.write('Ta\n')
		fdata2.write(f'{count}\n')
		if sys.argv[2] in ["c", "C"]:
			fdata2.write("Selective dynamics\n")
			fdata2.write("Cartesian\n")
			for item in pos_perfect:
			#print("{:12.9f} {:12.9f} {:12.9f}".format(cart_atoms[i][0], cart_atoms[i][1], cart_atoms[i][2]))
				fdata2.write(
					"{:12.12f} {:12.12f} {:12.12f} T T T\n".format(
						item[0] / lat_par, item[1] / lat_par, item[2] / lat_par
					)
				)


	return pos_perfect, count, basis, lx, ly, lz

'''##################################################################################
#------------------ Creating screw dislocation for randomize HEA -------------------#
##################################################################################'''	

def creating_screw_perfect(rand_atoms, count, basis, lx, ly, lz):	

	print ("Scaling of the cell (NX, NY, NZ) >>>", Nx, Ny, Nz)
	print ("Number of atoms >>>", count)
	print('>>>')
# A periodic array is quadrupole, if the vector d linking the two disloca-
# tions of opposite signs is equal to 1/2 (u1 +u2), where u1 and u2 are the periodicity
# vectors of the simulation cell. This ensures that every dislocation is a sym-
# metry center of the array: fixing, as a convention, the origin at a dislocation center,
# if a dislocation b is located at the position r, there will also be a dislocation b in −r.
#***************Position of dislocation line at Ci(X, Y, Z) for +/-b*********
	delta = 0.001     # to avoid on top of atom
	C1X = 0.25*lx + delta
	C2X = 0.75*lx + delta
	C1Y = 0.25*ly
	C2Y = 0.25*ly
	C3X = 0.75*lx + delta
	C4X = 0.25*lx + delta
	C3Y = 0.75*ly
	C4Y = 0.75*ly
	C1Z = C2Z = C3Z = C4Z = lz
	print('Generating Screw Dislocations at dislocation lines >>>')
	print("C1(X, Y, Z) for +b --> {:6.5f} {:6.5f} {:6.5f}".format(C1X, C1Y, C1Z) )
	print("C2(X, Y, Z) for -b --> {:6.5f} {:6.5f} {:6.5f}".format(C2X, C2Y, C2Z) )
	print("C3(X, Y, Z) for +b --> {:6.5f} {:6.5f} {:6.5f}".format(C3X, C3Y, C3Z) )
	print("C4(X, Y, Z) for -b --> {:6.5f} {:6.5f} {:6.5f}".format(C4X, C4Y, C4Z) )
	print('>>>')
	d1 = np.sqrt( (C1X-C2X)**2 + (C1Y-C2Y)**2 )
	d2 = np.sqrt( (C1X-C3X)**2 + (C1Y-C3Y)**2 )
	print("Dis b/w Dipoles d1, d2 --> {:6.5f}, {:6.5f}".format( d1, d2 ))
	jx = np.sqrt( np.sum(basis[0][:]*Nx)**2 )
	jy = np.sqrt( np.sum(basis[1][:]*Ny)**2 )
	print("length of supercell along jx, jy  --> {:6.5f}, {:6.5f}".format( jx, jy ))
	print("difference along jx-d1, (jx-d1)/2 --> {:6.5f}, {:6.5f}".format( jx-d1, (jx-d1)/2.0 ))
	print("difference along jy-d1, (jy-d2)/2 --> {:6.5f}, {:6.5f}".format( jy-d2, (jy-d2)/2.0 ))
#**************************************************************************

	burgers = lz # b = a * <111>/2	
	for l in range(count):
#--- (b/2pi)*tan**(-1)(y/x) ""DEF:: numpy.arctan2(x1, x2) x1 = y; x2 = x ""
		xx = rand_atoms[l][0]
		yy = rand_atoms[l][1]
		zz = rand_atoms[l][2]
#--- angle from X -> Y This is right
		theta1 = np.arctan2( (yy-C1Y ), (xx-C1X) ) 
		theta2 = np.arctan2( (yy-C2Y ), (xx-C2X) )
		theta3 = np.arctan2( (yy-C3Y ), (xx-C3X) )
		theta4 = np.arctan2( (yy-C4Y ), (xx-C4X) )
#--- Screw dislocation displacement in Z direction
		rand_atoms[l][2] = rand_atoms[l][2] + (burgers/(2.0*pivalue))*(theta1)
		rand_atoms[l][2] = rand_atoms[l][2] - (burgers/(2.0*pivalue))*(theta2)
		rand_atoms[l][2] = rand_atoms[l][2] + (burgers/(2.0*pivalue))*(theta3)
		rand_atoms[l][2] = rand_atoms[l][2] - (burgers/(2.0*pivalue))*(theta4)
#------------------------------ Writing to a POSCAR_Screw_unrelaxed file ------------------------------
	with open("POSCARrand_Screw_unrelaxed","w") as fdata1:
		fdata1.write(f'#{sys.argv[1]}\n')
		fdata1.write('{:6.6f}\n'.format(lat_par))
		fdata1.write("{:12.9f} {:12.9f} {:12.9f}\n".format(basis[0][0]*Nx/lat_par,basis[0][1]*Nx/lat_par,basis[0][2]*Nx/lat_par ))
		fdata1.write("{:12.9f} {:12.9f} {:12.9f}\n".format(basis[1][0]*Ny/lat_par,basis[1][1]*Ny/lat_par,basis[1][2]*Ny/lat_par ))
		fdata1.write("{:12.9f} {:12.9f} {:12.9f}\n".format(basis[2][0]*Nz/lat_par,basis[2][1]*Nz/lat_par,basis[2][2]*Nz/lat_par ))
		fdata1.write('Ta Hf Zr Nb Ti\n')
		fdata1.write('{0} {0} {0} {0} {0}\n'.format((int(count/5))) )
		if sys.argv[2] in ["c", "C"]:
			fdata1.write("Selective dynamics\n")			
			fdata1.write("Cartesian\n")	
			for i in range(len(rand_atoms) ):
			#print("{:12.9f} {:12.9f} {:12.9f}".format(cart_atoms[i][0], cart_atoms[i][1], cart_atoms[i][2]))
				fdata1.write("{:15.12f} {:15.12f} {:15.12f} T T T\n".format(rand_atoms[i][0]/lat_par, rand_atoms[i][1]/lat_par, rand_atoms[i][2]/lat_par ))
	
##################################################################################
#------------------------- Randomly distribute atoms for HEA---------------------#
##################################################################################	

def randomise_BCChea(rand_atoms, count, lat_vec):	
	rand_pos = []
	with open("POSCAR_perfect_rand","w") as fdata:
		for _ in range(int(1e3)):
			randomArrx = random.sample(range(count), count)

		if (count%5 == 0): 	
			fdata.write(f'#{sys.argv[1]}\n')
			fdata.write('{:6.6f}\n'.format(1.0))
			fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lat_vec[0][0]*Nx,lat_vec[0][1]*Nx,lat_vec[0][2]*Nx ))
			fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lat_vec[1][0]*Ny,lat_vec[1][1]*Ny,lat_vec[1][2]*Ny ))
			fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lat_vec[2][0]*Nz,lat_vec[2][1]*Nz,lat_vec[2][2]*Nz ))
			fdata.write('Ta Hf Zr Nb Ti\n')
			fdata.write('{0} {0} {0} {0} {0}\n'.format((int(count/5))) )
#------------------------- In Cartesian UNITS -------------------------#
			if sys.argv[2] in ["c", "C"]:
				fdata.write("Cartesian\n")
				print(
					f'# of atoms per element -> {count / 5}. File generated in Cartesian coordinates'
				)

				for i in range(len(rand_atoms) ):
					#print("{:12.9f} {:12.9f} {:12.9f}".format(rand_atoms[randomArrx[i]][0],rand_atoms[randomArry[i]][1],rand_atoms[randomArrz[i]][2] ) )
					fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(rand_atoms[randomArrx[i]][0],rand_atoms[randomArrx[i]][1],rand_atoms[randomArrx[i]][2] ) )
					rand_pos.append( rand_atoms[randomArrx[i]][:] )
#------------------------- In fractional/reduced UNITS ----------------#
			u = np.cross(lat_vec[1]*Ny, lat_vec[2]*Nz)
			v = np.cross(lat_vec[0]*Nx, lat_vec[2]*Nz)
			w = np.cross(lat_vec[0]*Nx, lat_vec[1]*Ny)
			V = np.array([ lat_vec[0]*Nx,lat_vec[1]*Ny,lat_vec[2]*Nz ] )
			print ("Volume of the cell::", LA.det(V) )
			Vx = np.inner(lat_vec[0]*Nx,u)
			Vy = np.inner(lat_vec[1]*Ny,v)
			Vz = np.inner(lat_vec[2]*Nz,w)			

			if sys.argv[2] in ["d", "D"]:
				fdata.write("Direct\n")
				print(
					f'# of atoms per element -> {count / 5}. File generated in Direct coordinates'
				)

				for i in range(len(rand_atoms) ):
					fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(np.dot(rand_atoms[randomArrx[i]],u)/Vx, np.dot(rand_atoms[randomArrx[i]],v)/Vy, np.dot(rand_atoms[randomArrx[i]],w)/Vz ) )
		else: 
			print(
				f'{count / 5} is not an integer number for equal composition -> HEAs not generated'
			)

		return rand_pos	

'''##################################################################################
#------------------------- WRITE POSCAR FILE ---------------------#
##################################################################################'''

def print_POSCAR(cart_atoms, count, lat_vec, NX, NY, NZ):	
	with open("POSCAR","w") as fdata:
		fdata.write(f'#{sys.argv[1]}\n')
		fdata.write('{:6.6f}\n'.format(1.0))
		fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lat_vec[0][0]*NX,lat_vec[0][1]*NX,lat_vec[0][2]*NX ))
		fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lat_vec[1][0]*NY,lat_vec[1][1]*NY,lat_vec[1][2]*NY ))
		fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(lat_vec[2][0]*NZ,lat_vec[2][1]*NZ,lat_vec[2][2]*NZ ))
		print("{:12.9f} {:12.9f} {:12.9f}".format(lat_vec[0][0]*NX,lat_vec[0][1]*NX,lat_vec[0][2]*NX ))
		print("{:12.9f} {:12.9f} {:12.9f}".format(lat_vec[1][0]*NY,lat_vec[1][1]*NY,lat_vec[1][2]*NY ))
		print("{:12.9f} {:12.9f} {:12.9f}".format(lat_vec[2][0]*NZ,lat_vec[2][1]*NZ,lat_vec[2][2]*NZ ))
		fdata.write('Ta\n')
		fdata.write(f'{count}\n')
#------------------------- In Cartesian UNITS -------------------------#
		if sys.argv[2] in ["c", "C"]:
			fdata.write("Cartesian\n")	
			for i in range(len(cart_atoms) ):
			#print("{:12.9f} {:12.9f} {:12.9f}".format(cart_atoms[i][0], cart_atoms[i][1], cart_atoms[i][2]))
				fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(cart_atoms[i][0], cart_atoms[i][1], cart_atoms[i][2]))
#------------------------- In fractional/reduced UNITS -------------------------#		
		u = np.cross(lat_vec[1]*NY, lat_vec[2]*NZ)
		v = np.cross(lat_vec[0]*NX, lat_vec[2]*NZ)
		w = np.cross(lat_vec[0]*NX, lat_vec[1]*NY)
		V = np.array([ lat_vec[0]*NX,lat_vec[1]*NY,lat_vec[2]*NZ ] )
		print ("Volume of the cell::", LA.det(V) )
		Vx = np.inner(lat_vec[0]*NX,u)
		Vy = np.inner(lat_vec[1]*NY,v)
		Vz = np.inner(lat_vec[2]*NZ,w)

		if sys.argv[2] in ["d", "D"]:
			fdata.write("Direct\n")		
			for i in range(len(cart_atoms) ):
				fdata.write("{:12.9f} {:12.9f} {:12.9f}\n".format(np.dot(cart_atoms[i],u)/Vx, np.dot(cart_atoms[i],v)/Vy, np.dot(cart_atoms[i],w)/Vz ) )	
 
			
def help():
	print('|-> Script to generate FCC or BCC supercell for HEAs.')
	print('|-> To execute just run python3 <bcc/fcc> <c/d> <sup/screw>.')
	print('|-> HEAs consists of five or more elements. The elements has already been typed into the')
	print('|-> script just change according to your needs, also lattice vectors should be')
	print('|-> equal and atomic composition should corresponds to integer multiple of atoms.')

if __name__ == '__main__':
	print (help())
	if (sys.argv[1] == 'bcc' and sys.argv[2] == 'c' and sys.argv[3] == 'screw'):
		pos_perfect, count, basis, lx, ly, lz = BCCgen_111_perfect()
		rand_pos = randomise_BCChea(pos_perfect, count, basis)			
		creating_screw_perfect(rand_pos, count, basis, lx, ly, lz)

	elif (
		sys.argv[1] in ['bcc', 'fcc']
		and sys.argv[2] in ['c', 'd']
		and sys.argv[3] == 'sup'
	):
		cart_atoms, count, lat_vec, SX, SY, SZ = HEAs_supercell()
		print_POSCAR(cart_atoms, count, lat_vec, SX, SY, SZ)
	else:
		print ("NOT VALID OPTION")
			
			
			
			
			
			
			
