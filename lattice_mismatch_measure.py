#!/usr/bin/env python3.6

'''

AUTHOR: ASIF IQBAL BHATTI
USAGE: python3 sys.argv[0]
DATED: 18/01/2020

'''

import numpy as np
import math, scipy, os, sys
from numpy import linalg as LA

np.set_printoptions(precision=6)

def local_lattice_distortion(a1,b1,c1):
	#print ("The lattice distortion in paracrystals is measured by the lattice distortion parameter g")
	#print (Back.YELLOW + "Wang, S. Atomic structure modeling of multi-principal-element alloys by the principle")
	#print (Back.YELLOW + "of maximum entropy. Entropy 15, 5536–5548 (2013).")
	#print ("")
	a=np.linalg.norm(a1); b=np.linalg.norm(b1); c=np.linalg.norm(c1)
	d = np.array([a,b,c])
	d_mean = np.mean(d); d_std = np.std(d)
	d_square_mean = (a**2 + b**2 + c**2)/3
	g = np.sqrt( d_square_mean/(d_mean)**2 - 1 )
	return g
	###
# Song, H. et al. Local lattice distortion in high-entropy alloys. Phys. Rev. Mater. 1, 23404 (2017).
# Senkov, O. N. & Miracle, D. B. Effect of the atomic size distribution on glass forming ability of amorphous metallic alloys. Mater. Res. Bull. 36, 2183–2198 (2001).
# Takeuchi, A. et al. Entropies in alloy design for high-entropy and bulk glassy alloys. Entropy 15, 3810–3821 (2013).	

def local_lattice_distortion_DEF1():
	#print ("The lattice distortion in paracrystals is measured by the lattice distortion parameter g")
	#print (Back.YELLOW + "Wang, S. Atomic structure modeling of multi-principal-element alloys by the principle")
	#print (Back.YELLOW + "of maximum entropy. Entropy 15, 5536–5548 (2013).")
	print ("+"*40,"HUME ROTHERY RULE","+"*40)
	C_i=C=0.2 ; r_avg = 0.0; del_sum=0.0
	elements = ["Nb", "Hf", "Ta", "Ti", "Zr"]
	eta = {
	"Nb" : 1.98,
	"Hf" : 2.08,
	"Ta" : 2.00,
	"Ti" : 1.76,
	"Zr" : 2.06, }
	
	print ("                      {element: atomic radius}")
	print (eta)
	
	for i in elements: 
		r_avg = r_avg + C * eta[i] 
	
	for j in elements:
		del_sum = del_sum + C * ( 1 - float(eta[j]) / r_avg )**2
	del_sum = 100 * np.sqrt(del_sum) 	
	print("HEA_atomic_size_mismatch: \u03B4={}".format(del_sum))
###
	
def local_lattice_distortion_DEF2():
	print ("Song, H. et al. Local lattice distortion in high-entropy alloys.")
	print ("Phys. Rev. Mater. 1, 23404 (2017).")
	print ("_____| Different definition of the atomic radius for the description ")
	print ("       of the local lattice distortion in HEAs")
	
	if not os.path.exists('POSCAR' and 'CONTCAR'):
		print (' ERROR: POSCAR & CONTCAR does not exist')
		sys.exit(0)
	print('Reading POSCAR and CONTCAR ... \n')
	
	x = []; y =[]; z=[]
	xp =[]; yp = []; zp = []; temp=0
	
	f = open('POSCAR','r')
	lines_poscar = f.readlines()
	f.close()
	
	f = open('CONTCAR','r')
	lines_contcar = f.readlines()
	f.close()
	
	sum_atoms = lines_poscar[6].split()  ### reading 7th lines for reading # of atoms
	sum_atoms = [int(i) for i in sum_atoms]
	sum_atoms = sum(sum_atoms)
	
	for i in lines_poscar:
		if "Direct" in i:
			lp=lines_poscar.index(i)
	for j in lines_contcar:
		if "Direct" in j:
			lc=lines_contcar.index(j)
			
	for i in range(sum_atoms):
		x, y, z    = lines_poscar[lp+1+i].split()
		xp, yp, zp = lines_contcar[lp+1+i].split()
		x = float(x); y = float(y); z = float(z)
		xp = float(xp); yp = float(yp); zp = float(zp)
		temp = temp + np.sqrt( (x-xp)**2 + (y-yp)**2 + (z-zp)**2 )
	temp = temp/sum_atoms
	print("local lattice distortion: \u0394d={}".format(temp))		

if __name__ == "__main__":	

	print ("+"*45,"HEA","+"*45)
	print ("-"*100)
	local_lattice_distortion_DEF1()
	print ("-"*100)
	local_lattice_distortion_DEF2()
	
	
	


