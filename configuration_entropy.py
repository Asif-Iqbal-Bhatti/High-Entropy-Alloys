#!/usr/bin/env python3

############################################################################
# The configuration entropy of a system with non-zero SRO can
# be estimated using the cluster variation method (CVM) in the 
# pair approximation [Tamm et al,]
############################################################################
import sys, os
import numpy as np	

k = 8.617333262145E-5 # Boltzmann constant eV/K
T = 300 # Temperature in Kelvin
z = 14 # Nearest neighbor 1st=4N, 2nd=3N, N is the number of lattice points
ci = [1/5,1/5,1/5,1/5,1/5] # concentration of each element in a alloy
print ("This script calcutes the ")
# Dictionary is defined for pairing probabilities. These values are taken from SRO.
dict = {'Nb-Ta':-0.05,
				'Nb-Ti':0.0,
				'Nb-Zr':0.1,
				'Nb-Hf':0.05,
				'Ti-Zr':0.125,
				'Ti-Ta':-0.025,
				'Ti-Hf':-0.15,
				'Ta-Zr':0.025,
				'Hf-Zr':-0.05,
				'Hf-Ta':0.1,
				'Nb-Nb':0.45,
				'Ta-Ta':0.475,
				'Ti-Ti':0.525,
				'Hf-Hf':0.525,
				'Zr-Zr':0.4}

for i in range (5):
	e1 = ( ci[i]*np.log(ci[i]) - ci[i] )
l1 = (z-1) * e1	
#print ( l1 )
for x in dict:
	e2 = ( dict[x]*np.log(dict[x]) - dict[x] )
l2 = z * e2
#print ( l2 )

# Final entropy of the system :
for x in dict:
	print ("{:6s} {:9.5f}".format( x, dict[x] ) )
print (  "-TS = -{:12.9f} meV for T=300K".format( (l1 - l2 + ( z/2 - 1 ))*k*1000*T    ))



