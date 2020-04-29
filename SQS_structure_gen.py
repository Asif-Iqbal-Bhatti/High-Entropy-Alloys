#!/usr/bin/python
# --------------------
# USAGE    : This script read foo.txt inputfile. The odering is important!!!
# NOTE     : A program to generate the SQS structure (either FCC or BCC) by using ATAT mcsqs code
# INPUTFILE: foo.txt is read by this code. foo.txt is the input parameters file.
#
# MODIFICATION of NANOHUB code by Asif Iqbal
# Created on: 30/04/2020 
#---------------------

import numpy as np
import sys
import os
import string
import subprocess

################# Reading from a file #########################
fo = open('foo.txt', 'r')
print "Name of the file: ", fo.name
print " "
elemcomp = fo.readline()
print "elemcomp: %s" % (elemcomp)
nat = fo.readline()
print "# of atoms: %d" % int((nat))
strtype = fo.readline()
print "Structure type is (BCC/FCC): %s" % str(strtype)
latparam = fo.readline()
print "Lattice parameter is: %f " % float((latparam))
cutoff = fo.readline()
print "Cutoff value are: %s" % (cutoff)
fo.close()
strtype = 'BCC'
################################################################

#################################################################

natoms=int(nat)
#print type(latparam)
value_ec=elemcomp.split()
elemtype=len(value_ec)/2
cutoff=cutoff.split()
#print cutoff[1]
element = []
composition = []
varstring = []

for i1 in range(2*elemtype):
   if (i1 % 2 ==0): 
       element.append(value_ec[i1]) 
   else:
       composition.append(value_ec[i1])
      # print composition
       
print element, composition
varlist = element
#varlist = list(string.ascii_uppercase)
#print varlist

ofile=open('input.in','w')
ofile.write('1  1  1  90   90  90  \n1  0  0  \n0  1  0  \n0  0  1 \n')

if (strtype == 'FCC'):
    for i in range(elemtype):
        varstring.append(varlist[i]+"="+composition[i])
    ofile.write('0.00   0.00   0.00 '+ ', '.join(varstring) + '\n')
    ofile.write('0.50   0.50   0.00 '+ ', '.join(varstring) + '\n')
    ofile.write('0.50   0.00   0.50 '+ ', '.join(varstring) + '\n')
    ofile.write('0.00   0.50   0.50 '+ ', '.join(varstring) + '\n')

if (strtype == 'BCC'):
    for i in range(elemtype):
        varstring.append(varlist[i]+"="+composition[i])

    ofile.write('0.00   0.00   0.00 '+ ', '.join(varstring) + '\n')
    ofile.write('0.50   0.50   0.50 '+ ', '.join(varstring) + '\n')
	
ofile.close()
latparam=float(latparam)

#os.system("corrdump -l=input.in -ro -noe -nop -clus -2=1 -3=1 -4=1")
os.system("corrdump -l=input.in -ro -noe -nop -clus -2='{}' -3='{}' -4='{}'".format(float(cutoff[0])/latparam, float(cutoff[1])/latparam, float(cutoff[2])/latparam)) 
os.system("mcsqs -n='{}' -l=input.in".format(natoms))
os.system("sort -k4 bestsqs-4.out > bestsqs-sorted.out")

bsqs = open("bestsqs-sorted.out","r")
bsqs1 = bsqs.readlines()
os.system("tail -{} bestsqs-sorted.out > tmp".format(natoms))
##############################################################
tmp4 = []
l = 0
tmp1 = open("tmp","r")
for line in tmp1:
    tmp3 = line.split()
    l += 1
    if (l % 10 == 0):
        tmp4.append(tmp3[3])
#print tmp4
###############################################################
pos = open("POSCAR","w")
ln=0 
for i in bsqs1:
    if(ln==2):
       pos.writelines(element)
       pos.write("\n1.0\n")
    if(ln>2):
       j=i.split()
       posstring=[str(float(j[k])*latparam)+' ' for k in range(3)]
       for k in range(elemtype):
           try:
              if(j[3] == varlist[k]):
                 #type = j[3]
                 posstring.append(varlist[k])
                 #print type
           except:
               pass
       pos.writelines(posstring)
       pos.write("\n")
       if(ln==5):  
          pos.writelines([elem+" " for elem in tmp4])
          pos.write("\n")
          pos.writelines([str(float(comp)*natoms)+" " for comp in composition])
          pos.write("\nCartesian\n")
    ln=ln+1   
pos.close()
os.system("column -t POSCAR > POSCAR-1 && mv POSCAR-1 POSCAR")

# for testing
#print 'formula = %s' % elemcomp
#print 'Total Atoms = %d' % natoms
#print 'Str Type = %s' % strtype
#print 'Lattice = %s' % latparam
#print 'Nr of Element Type = %d' % elemtype

sys.exit()

