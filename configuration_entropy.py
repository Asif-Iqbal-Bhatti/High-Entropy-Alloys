#!/usr/bin/env python3

############################################################################
# The configuration entropy of a system with non-zero SRO can
# be estimated using the cluster variation method (CVM) in the 
# pair approximation [Tamm et al,]
# NOT SURE ABOUT THIS METHOD. 
############################################################################
import sys, os
import numpy as np	

k = 8.617333262145E-5 # Boltzmann constant eV/K
T = 300 # Temperature in Kelvin
z = 14 # Nearest neighbor 1st=4N, 2nd=3N, N is the number of lattice points
ci = [1/5,1/5,1/5,1/5,1/5] # concentration of each element in a alloy

# Dictionary is defined for pairing probabilities. These values are taken from SRO.
e1 = 0
e2 = 0
dd = []
os.remove("sro_tmp.dat")
dict = {'Nb-Ta':0.05,
				'Nb-Ti':0.0002,
				'Nb-Zr':0.1,
				'Nb-Hf':0.05,
				'Ti-Zr':0.125,
				'Ti-Ta':0.025,
				'Ti-Hf':0.15,
				'Ta-Zr':0.025,
				'Hf-Zr':0.05,
				'Hf-Ta':0.1,
				'Nb-Nb':0.45,
				'Ta-Ta':0.475,
				'Ti-Ti':0.525,
				'Hf-Hf':0.525,
				'Zr-Zr':0.4}
##################################################################################				

sro_tmp = os.popen("sqsgenerator alpha sqs POSCAR --weight=1,0.5 --verbosity=3 " ).read()
with open("sro_tmp.dat", 'w') as output:
	for line in sro_tmp:
		output.write(line)
with open("sro_tmp.dat") as fp:
	line = fp.readline()
	cnt = 1
	while line:
		#print("{}".format(line.strip(" ").strip("\n").strip(";").split("=") ))
		dd.append( line.strip(" ").strip("\n").strip(";").split("=") )
		cnt +=1
		line = fp.readline()
##################################################################################

for item in ci:
	e1 = e1 + (item * np.log(item) - item)
l1 = (z-1) * e1	

#for x in dict:
#	e2 = e2 + ( dict[x]*np.log(dict[x]) - dict[x] )
#l2 = z/2 * e2

for i in range(cnt-2):
	e2 = e2 + ( float(dd[i+1][1]) * np.log(float(dd[i+1][1])) - float(dd[i+1][1]) )
	#print (dd[i+1][0], float(dd[i+1][1]))
l2 = z/2 * e2	

# Final entropy of the system :
for x, value in dict.items():
	print("{:6s} {:9.5f}".format(x, value), end='\r')
print (  "T*S = {:12.5f} meV at T=300K".format( (l1 - l2 + ( z/2 - 1 ))*k*1000*T    ))

