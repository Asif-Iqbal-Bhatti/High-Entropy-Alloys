#!/usr/bin/env python3

'''
##############################################################################
# ML error analysis
# AIB_EM for VASP.6.4.1
##############################################################################
'''

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

os.system("grep T OSZICAR > frameVsEne.dat")
file_path = 'frameVsEne.dat'
df = pd.read_csv(file_path, delim_whitespace=True, header=None)
#df.index = range(1, len(df) + 1)
print(df)
df.to_csv('frameVsEne.csv')
min_config = df[4].idxmin()
print('MINIMUM INDEX :: ', min_config,'\n', df.iloc[min_config].T)

#######################  BEEF #######################
os.system("cat ML_LOGFILE | grep BEEF > BEEF.dat")
with open('BEEF.dat', 'r') as file:
    header_line = file.readlines()[11].strip()
columns1 = header_line.split()[2:]
#print(columns1)
data = np.genfromtxt('BEEF.dat', skip_header=14, usecols=range(1, 8),)
df1 = pd.DataFrame(data, columns=columns1)
df_filtered = df1[df1['nstep'] >= 100]
df_filtered.plot(x='nstep')
plt.savefig('BEEF.pdf', bbox_inches='tight', dpi=300)
plt.clf()
    
####################### ERR #######################
os.system("cat ML_LOGFILE | grep ERR > ERR.dat")
with open('ERR.dat', 'r') as file:
    header_line = file.readlines()[8].strip()
columns2 = header_line.split()[2:]
#print(columns2)
data = np.genfromtxt('ERR.dat', skip_header=11, usecols=range(1, 5),)
df2 = pd.DataFrame(data, columns=columns2)
df_filtered = df2[df2['nstep'] >= 100]
df_filtered.plot(x='nstep')
plt.savefig('ERR.pdf', bbox_inches='tight', dpi=300)
plt.clf()


