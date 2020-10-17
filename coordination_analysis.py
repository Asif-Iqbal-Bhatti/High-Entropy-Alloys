#!/usr/bin/env python3
############################################################################
# USAGE  :: python3 coordination_analysis.py 
# Author :: Asif Iqbal
# DATED  :: 17/10/2020
# NB		 :: POSCAR should be in Cartesian coordinates.
# Calculate the coordination around the dislocation line and
# counts the number and types of atoms within the cutoff radius.
# This code reads the POSCAR file with atomic types appended to the last line.
############################################################################

import numpy as np
from scipy import stats 
import os, sys, random, subprocess, shutil
from matplotlib import pyplot as plt
import pandas
from collections import Counter

def read_poscar():
	pos = []; kk = []; lattice = []; atom_index = []; sum = 0; dict = {}
	file = open('POSCAR','r')
	firstline  = file.readline() # IGNORE first line comment
	alat       = float( file.readline() )# scale
	Latvec1    = file.readline().split(); #print("{:9.6f} {:9.6f} {:9.6f}".format(float(Latvec1[0]),float(Latvec1[1]),float(Latvec1[2])))
	Latvec2    = file.readline().split(); #print("{:9.6f} {:9.6f} {:9.6f}".format(float(Latvec2[0]),float(Latvec2[1]),float(Latvec2[2])))
	Latvec3    = file.readline().split(); #print("{:9.6f} {:9.6f} {:9.6f}".format(float(Latvec3[0]),float(Latvec3[1]),float(Latvec3[2]))) 
	elementtype= file.readline(); #print ("{}".format(elementtype.split() ))
	atomtypes  = file.readline(); #print ("{}".format(atomtypes.split() ))
	Coordtype  = file.readline().split()
	if (Coordtype[0] == 'Direct' or Coordtype[0] == 'direct'): exit("Convert to Cartesian, first!")
	nat = [int(i) for i in atomtypes.split()]
	for i in nat: sum = sum + i; n_atoms = sum				
	for x in range(int(n_atoms)):
		coord = file.readline().rstrip(" ").split()
		coord = [ str(i) for i in coord[0:4] ]
		pos 	= pos + [coord]
	file.close() 
	return n_atoms,pos,firstline,alat,Latvec1,Latvec2,Latvec3,elementtype,atomtypes,Coordtype

def replicate_cell(pos,n_atoms,Latvec1,Latvec2,Latvec3):
	Nx,Ny,Nz = 2,2,2; mag_atoms_pos = []; mm = []; n = []
	Latvec1  = [ float(Latvec1[0]), float(Latvec1[1]), float(Latvec1[2]) ]
	Latvec2  = [ float(Latvec2[0]), float(Latvec2[1]), float(Latvec2[2]) ]
	Latvec3  = [ float(Latvec3[0]), float(Latvec3[1]), float(Latvec3[2]) ]
	Latvect  = np.array( [Latvec1[:], Latvec2[:], Latvec3[:]] )
	for l in range(n_atoms):
		ll = [float(pos[l][0]), float(pos[l][1]), float(pos[l][2]), pos[l][3] ]
		mm.append( ll )	
	for l in mm:
		for i in range(0,Nx,1):            
			for j in range(0,Ny,1):
				for k in range(0,Nz,1):          
					cartesian_basis = np.inner(Latvect.T, np.array( [i,j,k] ) )
					mag_atoms_pos.append( cartesian_basis + l[0:3] ) 	
					n.append( l[3:4] )
	Nx, Ny, Nz = 1,1,1		
	ff = open("POSCAR_222", 'w')
	ff.write("POSCAR_{}x{}x{}\n".format(Nx,Ny,Nz))
	ff.write("1.0\n")
	ff.write("{:9.7f} {:9.7f} {:9.7f}\n".format(Nx*Latvec1[0], Ny*Latvec1[1], Nz*Latvec1[2] ) )
	ff.write("{:9.7f} {:9.7f} {:9.7f}\n".format(Nx*Latvec2[0], Ny*Latvec2[1], Nz*Latvec2[2] ) )
	ff.write("{:9.7f} {:9.7f} {:9.7f}\n".format(Nx*Latvec3[0], Ny*Latvec3[1], Nz*Latvec3[2] ) )
	ff.write("Hf   Nb   Ta   Ti   Zr\n")
	ff.write(5*"{:6d}".format( int(len(mag_atoms_pos)/5) ))
	ff.write("\nCartesian\n")
	for g in range(len(mag_atoms_pos)):
		ff.write( "{:12.7f} {:12.7f} {:12.7f}\n" .format(mag_atoms_pos[g][0], mag_atoms_pos[g][1], mag_atoms_pos[g][2] ))
	ff.close()
	return  mag_atoms_pos, n
	
def coordination_analysis(n_atoms, pos):
	# First select the atom around which to measure the coordination.
	# I've used (x-x0)**2 + (y-y0)**2 + (z-z0)**2 < R**2; (r-r0)**2 < R**2
	r0 = 313;	cnt = 0; rcutoff = 3.0; bar_graph = [];	
	for x in range(0, int(pos), 1):
		i = float(pos[x][0]) - float(pos[r0][0]) 
		j = float(pos[x][1]) - float(pos[r0][1]) 
		k = float(pos[x][2]) - float(pos[r0][2]) 
		if ( ( i*i + j*j + k*k ) <= rcutoff**2 ): #FILTER
			cnt +=1; bar_graph.append( pos[x][3] )	
			print ( "{:4d} {}".format(x,pos[x][:]) )
	print ("Coordination # {}".format( cnt) )
	element_counts = Counter(bar_graph)
	return element_counts, bar_graph

# Minimum image convention needs to be implemented, to fully account the
# coordination. It does not work by replicating the system.

def coordination_analysis_replicate_cell(n_atoms, pos, n):
	# First select the atom around which to measure the coordination.
	# I've used (x-x0)**2 + (y-y0)**2 + (z-z0)**2 < R**2; (r-r0)**2 < R**2
	r0 = 868;	cnt = 0; rcutoff = 3.0; bar_graph = [];	
	for x in range(0, len(pos), 1):
		i = float(pos[x][0]) - float(pos[r0][0]) 
		j = float(pos[x][1]) - float(pos[r0][1]) 
		k = float(pos[x][2]) - float(pos[r0][2]) 
		if ( ( i*i + j*j + k*k ) < rcutoff**2 ): #FILTER
			cnt +=1; print ( "{:4d} {} {}".format(x, pos[x][:],n[x] ))
	print ("Coordination # {}, atoms {}".format( cnt, len(pos)) )
	element_counts = Counter(bar_graph)
	return element_counts, bar_graph
	
def plot_bar_from_counter(counter, bar_graph, ax=None):
	if ax is None:
		fig = plt.figure()
		ax = fig.add_subplot(111)
	frequencies = counter.values()
	names = counter.keys()
	#ax.bar(names, frequencies, align='center')
	#ax.bar(names, frequencies)
	plt.hist(bar_graph, bins=20, facecolor='green', alpha=0.75)
	plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':300})
	plt.gca().set(title='Frequency Histogram', ylabel='Frequency');
	plt.show()
	return ax

#-------------------------------------MAIN ENGINE-----------------------------------------
if __name__ == "__main__":
	n_atoms,pos,firstline,alat,Latvec1,Latvec2,Latvec3,elementtype,atomtypes,Coordtype = read_poscar();
	Burgers = float(Latvec3[2]); # B = a/2<111> which is length along Z=<111>
	alat = (1 * Burgers ) / np.sqrt(3)
	
	print("SYSTEM detected={}, atoms={} alat={:5.5f}\n".format(elementtype.split(), n_atoms, alat), end='\t' )
	print("{:10s} {:20.12s} {:15.12s}".format("atom#", "position", "atom type"))
	
	mag_atoms_pos, n = replicate_cell(pos,n_atoms,Latvec1,Latvec2,Latvec3)	
	element_counts, bar_graph = coordination_analysis_replicate_cell(n_atoms, mag_atoms_pos, n)
	
	#element_counts = coordination_analysis(n_atoms, pos)
	#plot_bar_from_counter(element_counts, bar_graph)
	
	
