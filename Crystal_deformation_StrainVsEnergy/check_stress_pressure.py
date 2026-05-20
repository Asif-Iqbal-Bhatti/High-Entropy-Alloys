#!/usr/bin/env python3

'''
A simple script to extract stress tensor from the vasp xml file and 
compute the pressure and correction to the stress tensor.

This script uses ase module to exttract data instead of writing the code
from scratch.
https://wiki.fysik.dtu.dk/ase/_modules/ase/io/vasp.html#read_vasp_xml
https://wiki.fysik.dtu.dk/ase/_modules/ase/atoms.html#Atoms.get_forces

'''
#-------------------------------------------------------------------------------#
import numpy as np
import ase.units
#from ase import Atoms
import ase.io
from ase.units import GPa
from numpy import dot, diag, ones, reshape, linspace, array, mean
from math import acos, pi, cos, sin, sqrt
CRED = '\033[91m';CEND = '\033[0m'
CYEL = '\033[33m'; CEND = '\033[0m'
CPIN = '\033[46m';
#-------------------------------------------------------------------------------#
# TODO: Needs verification if this ase code is correct the units are
# 1 eV/Angstrom3
I = np.identity(3)
file = ase.io.read('vasprun.xml')

print ("\n{:*^80s}\n".format("Reading directly vasprun.xml by ase code"))
print ("Reading vasprun.xml file. Stress tensor is given in matrix notation:")
print ("Volume -> {:6.3f} A**3".format(file.get_volume() ))

print (CRED + "Length&Angles->{}".format(file.get_cell_lengths_and_angles() ) + CEND)
print (CRED + "Lattice vectors(A)::\n", end="" + CEND)
print ("{}".format(np.mat(file.get_cell()) ), end="\n")
print ("="*80)

# stress matrix is already multiplied by *-1 in VASP
stress = file.get_stress(voigt=False) # this can be switched
print (CYEL+ "Stress tensor in Matrix notation (eV/Angstrom3):\n" + CEND)
print (CYEL+ "{}\n".format(stress) + CEND)

p = -1*(stress.trace())/3 # Hydrostatic Pressure
hyd_stress = (stress.trace())/3 # Hydrostatic Stress

print ("pressure => {:5.3f} GPa".format(p*160.21766208))
print ("Stress-pressure =>\n",stress - dot(p,I) )


#-------------------------------------------------------------------------------#
'''
# This is the manual demonstration of how ase reads vasprun.xml file
# One thing to remember in vasp stress units are in kb which can be 
# converted to GPa by 0.1*GPa

  param :: pressure = - stress
	1kbar = 0.1 * GPa
	In vasprun.xml stress units are given in kbar
'''
#-------------------------------------------------------------------------------#
print ("{:_^80s}".format("Example my own conversion"))
print ("In voigt notation:")
stress1 = np.zeros((3, 3), dtype=float)
print (GPa)
sblocks = [
   [ -0.05224891,   0.00659668,   0.00785320],
   [  0.00659642,  -0.07778282,   0.00069340],
   [  0.00785295,   0.00069324,  -0.07633123] ]
s = np.array(sblocks)	 
if s is not None:
	for i, v in enumerate(s):
		stress1[i] = np.array([float(val) for val in v[:]])
	stress1 = stress1 * (-0.1)	
	print(stress1)	
	
	stress1 = stress1.reshape(9)[[0, 4, 8, 5, 2, 1]]
	p = - mean(stress1[:3]) # Hydrostatic Pressure	
	
	print ("pressure => {:5.3f} GPa".format(p))
	print ("Stress - pressure =>\n", stress1 - p )
	
	
	
	
	
	
	
	
	
