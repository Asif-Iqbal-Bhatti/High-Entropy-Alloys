#!/usr/bin/env python3
	
'''
############################################################################
# USAGE  :: python3 sys.argv[0] 
# Author :: Asif Iqbal
# DATED  :: 14/12/2020
############################################################################
'''

import numpy as np
import os, sys, random, subprocess, shutil, math
from termcolor import colored
subscript = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
from sympy import *
import sympy as sp
sp.init_printing()
	
'''
#Composition is given by this relation.
#    (ABCD)_(1-x)/4 . A_(x), where x=[0, 1]
#
#The formation energy can be defined as:
#		 H_form = E_HEA/N - SUM_i {Ci * Eref_i/N}; i = {1 ... 5}
#
#    The elements and energies are in 1-1 correspondence
#    The bulk phase energies are reported for two atoms per unit cell. 
'''
E_HEA, N, C_i, Eref_i, C, n, i = sp.symbols('E_HEA, N, C_i, Eref_i, C, n, i')

H_form = E_HEA/N - sp.summation( C_i * Eref_i / N, (i, 0, n) );
pprint(H_form)

A, ABCD, x = sp.symbols('A_x, (ABCD)_(1-x)/4, x', integer=True)
comp_exp =  ABCD * A ;
print(comp_exp, "; x[0, 1]")
print("="*50)
HEA_elem   = ['A', 'B', 'C', 'D', 'E']
HEA_energy = [-25.475932, -20.433322, -23.629510, -17.041916, -15.673653]	

EA = []; c = 0

for x in np.arange(0, 1.1, 0.1):

	X = round ( (1-x)/4, 3 )
	Y = round ( x, 1 ) 
	
	EA = HEA_elem[0]+str(X)+HEA_elem[1]+str(X)+HEA_elem[2]+str(X)+\
	HEA_elem[3]+str(X)+HEA_elem[4]+str(Y)
	
	c = 4*X + Y	
	print(("{:30.35s} | Tot= {} |".format( EA.translate(subscript), c)), end=" ")
	print("Ci: [{:5.1f}%, {:5.1f}%]".format( X*100, Y*100), end=" ")
	
	JJ = HEA_energy[0]*X + HEA_energy[1]*X + HEA_energy[2]*X + HEA_energy[3]*X + HEA_energy[4]*Y
	
	print("EForm: {:8.3f}".format( -1267.01046136/125 - JJ/2) )
print("="*50)


