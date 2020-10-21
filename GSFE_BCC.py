#!/usr/bin/env python3
#---------------------------------------------
# Generate coordinate file for generalized stacking fault energy calculation for Ta_{110}<111>
# bcc burger vector = 1/2*[111] slip plane = (110): [111]:x, [112]:y, [110]:z ###
#---------------------------------------------
import os
import re
import sys
import math

# Lattice parameter of Ta bcc
a = 3.319

def get_bcc(a,t,nx=1,ny=1,nz=1):
  xa = []; ya = []; za = [];tz = [];lz=[]; step = 10
	#components of a vector
  ax = a*math.sqrt(3)/2.
  ay = a*math.sqrt(6)
  az = a*math.sqrt(2)
	#duplication of cell
  bx,by,bz = ax*nx,ay*ny,az*nz
  ux = (ax/step)*t; print("{:2d} {:6.4f}".format(t, ux) )
  layer=0
  for k in range(int(nz/2) ):
    for j in range(ny):
      for i in range(nx):
        xa.append((0/6.+i)*ax); ya.append((0/6.+j)*ay); za.append((0/6.+k)*az);lz.append(layer)
        xa.append((0/6.+i)*ax); ya.append((1/2.+j)*ay); za.append((1/2.+k)*az);lz.append(layer);layer+=1
        xa.append((6/9.+i)*ax); ya.append((1/3.+j)*ay); za.append((0/2.+k)*az);lz.append(layer)
        xa.append((6/9.+i)*ax); ya.append((5/6.+j)*ay); za.append((1/2.+k)*az);lz.append(layer);layer+=1
        xa.append((3/9.+i)*ax); ya.append((1/6.+j)*ay); za.append((1/2.+k)*az);lz.append(layer)
        xa.append((3/9.+i)*ax); ya.append((6/9.+j)*ay); za.append((0/2.+k)*az);lz.append(layer)
  for k in range(int(nz/2),nz):
    for j in range(ny):
      for i in range(nx):
        xa.append((0/6.+i)*ax+ux); ya.append((0/6.+j)*ay); za.append((0/6.+k)*az);lz.append(layer)
        xa.append((0/6.+i)*ax+ux); ya.append((1/2.+j)*ay); za.append((1/2.+k)*az);lz.append(layer);layer+=1
        xa.append((6/9.+i)*ax+ux); ya.append((1/3.+j)*ay); za.append((0/2.+k)*az);lz.append(layer)
        xa.append((6/9.+i)*ax+ux); ya.append((5/6.+j)*ay); za.append((1/2.+k)*az);lz.append(layer);layer+=1
        xa.append((3/9.+i)*ax+ux); ya.append((1/6.+j)*ay); za.append((1/2.+k)*az);lz.append(layer)
        xa.append((3/9.+i)*ax+ux); ya.append((6/9.+j)*ay); za.append((0/2.+k)*az);lz.append(layer)

# Wraping the atoms inside the unit cell. One fault per unit cell. 
  vacuum = 0.0 # It is better to first relax the structure and then introduce vacuum
  for i in range(len(xa)):
    xa[i] = ( (xa[i] + bx)/bx - int((xa[i] + bx)/bx) )* bx
    ya[i] = ( (ya[i] + by)/by - int((ya[i] + by)/by) )* by
    za[i] += vacuum/2.0
  bz += vacuum
  return xa,ya,za,bx,by,bz,lz, layer

def gen_poscar(xa,ya,za,box,tau):
  fout = open("POSCAR"+str(tau),"w")
  fout.write("Ta bcc (110): [111]:x, [112]:y, [110]:z\n")
  fout.write("1.0\n")
  fout.write("{:22.16f}  {:22.16f}  {:22.16f}\n".format(box[0],0,0))
  fout.write("{:22.16f}  {:22.16f}  {:22.16f}\n".format(0,box[1],0))
  fout.write("{:22.16f}  {:22.16f}  {:22.16f}\n".format(0,0,box[2]))
  fout.write("%d\n"%len(xa))
  fout.write("Selective Dynamics\n")
  fout.write("Cart\n")
  for i in range(len(xa)):
    fout.write("%22.16f %22.16f %22.16f F F T\n"%(xa[i],ya[i],za[i]))
  fout.close()
  return len(xa)
# ---------------------------------------------------------
usage="""
        Usage: bcc_111_gsfe_curve.py a nx ny nz       
               a  - equilibrium lattice parameter
               nx,ny,nz - periodicity in x,y and z
        Example: bcc_111_gsfe_curve.py 3.315 1 1 6 
"""
# ---------------------------------------------------------
if __name__ == '__main__':

	if len(sys.argv) > 4:
		a=float(sys.argv[1])
		nx=int(sys.argv[2])
		ny=int(sys.argv[3])
		nz=int(sys.argv[4])

#   step size should be the same as the range of tau
		for tau in range(11): 
			xa,ya,za,bx,by,bz,lz,layer = get_bcc(a,tau,nx,ny,nz)
			gen_poscar(xa,ya,za,[bx,by,bz],tau)
		print ("# of layers -> ", layer)			
		print ("--> Files are generated")			
	else:
		print ("Error: wrong number of arguments!!!")
		print (usage)
		
