#!/usr/bin/env python3

'''
###
# ML error analysis
# AIB_EM for VASP.6.4.1
###
'''

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


file_path = 'frameVsEne.dat'
df = pd.read_csv(file_path, delim_whitespace=True, header=None)
print(df)
print(df[4].idxmin())

######
os.system("cat ML_LOGFILE | grep BEEF > BEEF.dat")
with open('BEEF.dat', 'r') as file:
    header_line = file.readlines()[11].strip()
columns = header_line.split()[2:]
print(columns)

data = np.genfromtxt('BEEF.dat', skip_header=14, usecols=range(1, 8),)
df = pd.DataFrame(data,)
df.columns = columns
print(df)
df.plot(x='nstep')

######
os.system("cat ML_LOGFILE | grep ERR > ERR.dat")
with open('ERR.dat', 'r') as file:
    header_line = file.readlines()[8].strip()
columns = header_line.split()[2:]
print(columns)

data = np.genfromtxt('ERR.dat', skip_header=11, usecols=range(1, 5),)
df = pd.DataFrame(data,)
df.columns = columns
print(df)
df.plot(x='nstep')
plt.show()