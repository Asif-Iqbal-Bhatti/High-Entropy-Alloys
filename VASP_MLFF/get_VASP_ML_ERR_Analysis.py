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

os.system("grep T OSZICAR > frameVsEne.dat")
df = pd.read_csv(file_path, delim_whitespace=True, header=None)
minV = df[4].idxmin()
print(minV)
print((df.iloc[minV]).T)

######
os.system("cat ML_LOGFILE | grep BEEF > BEEF.dat")
with open('BEEF.dat', 'r') as file:
    header_line = file.readlines()[11].strip()
columns = header_line.split()[2:]

data = np.genfromtxt('BEEF.dat', skip_header=14, usecols=range(1, 8),)
df = pd.DataFrame(data,)
df.columns = columns
df.plot(x='nstep')
plt.savefig('BEEF.png', bbox_inches='tight', dpi=300)
plt.clf()

######
os.system("cat ML_LOGFILE | grep ERR > ERR.dat")
with open('ERR.dat', 'r') as file:
    header_line = file.readlines()[8].strip()
columns = header_line.split()[2:]

data = np.genfromtxt('ERR.dat', skip_header=11, usecols=range(1, 5),)
df = pd.DataFrame(data,)
df.columns = columns
df.plot(x='nstep')
plt.savefig('ERR.png', bbox_inches='tight', dpi=300)
plt.clf()

######
os.system("grep BEEF ML_LOGFILE|grep -v '#'|awk '{print $2, $6}' > CTIFOR.dat")
data = np.genfromtxt('CTIFOR.dat')
df = pd.DataFrame(data,)
df.plot(x=0)

