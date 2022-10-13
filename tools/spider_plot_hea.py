#!/usr/bin/env python3

'''
#===========================
AUTHOR:: ASIF IQBAL -> @AIB_EM
USAGE:: spider plot for HEA
#===========================
'''

import pandas as pd
import plotly.graph_objects as go
from matplotlib import pyplot as plt, patches
import plotly.offline as pyo
import numpy as np

data2nn = 'Ti$_{x=0.20}$_2NN'
data3nn = 'Ti$_{x=0.20}$_3NN'

#============= DATA_2NN
data_2NN = pd.read_csv('VASPRun_x0.20_800_1st_2nd.csv')
data_2NN = data_2NN.sort_values(by=['TE/at'], ascending=True)
slow2NN = data_2NN.loc[data_2NN['TE/at'].idxmin()]
shig2NN = data_2NN.loc[data_2NN['TE/at'].idxmax()]

at_pair = list(data_2NN.columns[1:-3])
at_pair = [*at_pair, at_pair[0]]

plow2NN = list(slow2NN[1:-3])
plow2NN = [*plow2NN, plow2NN[0]]

phig2NN = list(shig2NN[1:-3])
phig2NN = [*phig2NN, phig2NN[0]]

#========== DATA_3NN
data3NN = pd.read_csv('VASPRun_x0.20_800_1st_2nd_3rd.csv')
data3NN = data3NN.sort_values(by=['TE/at'], ascending=True)
slow3NN = data3NN.loc[data3NN['TE/at'].idxmin()]
shig3NN = data3NN.loc[data3NN['TE/at'].idxmax()]

plow3NN = list(slow3NN[1:-3])
plow3NN = [*plow3NN, plow3NN[0]]

phig3NN = list(shig3NN[1:-3])
phig3NN = [*phig3NN, phig3NN[0]]

#========== RADAR CHART
label_loc = np.linspace(start=0, stop=2 * np.pi, num=len(plow2NN))
plt.figure(figsize=(6, 6))
plt.subplot(polar=True)

plt.scatter( 0.0 , -0.4, s=35000 , linewidth=1.5, facecolors='none', edgecolors='red', linestyle='--') 
plt.plot(label_loc, plow2NN, label=f'{data2nn}_low', linewidth=1.2, linestyle='solid', marker='o')
plt.plot(label_loc, phig2NN, label=f'{data2nn}_high', linewidth=1.2, linestyle='solid', marker='o')

plt.plot(label_loc, plow3NN, label=f'{data3nn}_low', linewidth=1.2, linestyle='solid', marker='s')
plt.plot(label_loc, phig3NN, label=f'{data3nn}_high', linewidth=1.2, linestyle='solid', marker='s')

lines, labels = plt.thetagrids(np.degrees(label_loc), labels=at_pair)

plt.title('SRO comparison of low and high configurations Ti$_{x=0.20}$', size=14, y=1.05)
plt.text(2.6, 0.6, '(b)', fontsize = 22)
plt.legend(loc='upper right', bbox_to_anchor=(1.2, 1.0), prop={'size': 7})

plt.yticks(color="grey", size=8)
plt.xticks(label_loc, at_pair, color='black', size=14)
plt.ylim([-0.4,0.3])

plt.savefig('SRO_Tix0.20_plot.png', dpi=600,bbox_inches='tight')

