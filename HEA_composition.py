#!/usr/bin/env python3

############################################################################
# USAGE  :: python3 sys.argv[0] 
# Author :: Asif Iqbal
# DATED  :: 14/16/2022
############################################################################

import numpy as np
import os, sys, random, subprocess, shutil, math
from termcolor import colored
subscript = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

#>>> Composition is given by this relation.
#    (HfNbTaZr)_(1-x)/4 . Ti_(x), where x=[0, 1]
#
#>>> The formation energy can be defined as:
#		 H_form = E_HEA/N - SUM_i {Ci * Eref_i/N}; i = {1 ... 5}
#
#    The elements and energies are in 1-1 correspondence
#    The energies are reported for two atoms per unit cell. 

HEA_elem   = ['Hf', 'Nb', 'Ta', 'Zr', 'Ti']
HEA_energy = [-25.475932, -20.433322, -23.629510, -17.041916, -15.673653]	

EA = []; c = 0

for x in np.arange(0, 1.1, 0.1):

	A = round ( (1-x)/4, 3 )
	B = round ( x, 1 ) 
	
	EA=HEA_elem[0]+str(A)+HEA_elem[1]+str(A)+HEA_elem[2]+str(A)+HEA_elem[3]+str(A)+HEA_elem[4]+str(B)
	c = 4*A + B
	
	#shutil.rmtree(EA, ignore_errors=True) #overwrite a directory
	#os.mkdir(EA)
	
	print(colored("{:30.35s} > sum={}".format( EA.translate(subscript), c), 'red'), end=" ")
	print("Concentration Ci = {}%, {}%".format( A*100, B*100))
	
	JJ = HEA_energy[0]*A + HEA_energy[1]*A + HEA_energy[2]*A + HEA_energy[3]*A + HEA_energy[4]*B
	
	print("Formation energy is {:.6f}".format( -1267.01046136/125 - JJ/2) )



