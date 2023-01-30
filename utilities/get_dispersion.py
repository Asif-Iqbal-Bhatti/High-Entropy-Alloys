#!/usr/bin/env python3

# AUTHOS :: Asif
# USAGE  :: get dispersion 'Edisp' from OUTCAR file in a pythonic way
# DATED  :: 30/01/2023

import os,sys, re, subprocess, gzip
from monty.io import zopen

path = os.path.join('03_ISIF3snsb','OUTCAR.gz')
print(path)


edisps=[]
with zopen(path, 'rb') as f:
	for line in (f):
		line = line.decode().strip()
		match = re.search(r'Edisp.*?([-\d]+\.\d+)', line)
		if match:
			edisps.append(float(match.group(1)))
print(edisps[-1])
