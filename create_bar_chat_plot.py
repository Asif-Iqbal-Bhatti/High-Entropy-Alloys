#!/usr/bin/env python

'''
# NOTE:: skript voor bar plot X vs Y
# https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.plot.bar.html
# FIRST CONVERT => *.txt or *.dat into csv file
'''

import pandas as pd
import matplotlib.pyplot as plt
	
a = {}
subscript = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

#======  MAIN WORK
dg = pd.read_csv('X_vs_Y.csv')
naam = list(dg.columns)

for j in range(len(dg)):
	a[j] = [dg[naam[0]][j].translate(subscript)] + [dg[naam[1]][j]]

df = pd.DataFrame(a).T
df.columns = naam
df = df.sort_values(by=naam[1], ascending=False)

df[naam[1]] = df[naam[1]].astype(float)
print(df)

#====== WRITING TO A FILE
ax = df.plot.bar(color={naam[1]: "r"})
ax.set_xticklabels(df['Formula'], rotation=60, ha='right',fontweight='light',fontsize='medium')
ax.set_ylabel(naam[1], rotation=90,fontweight='light',fontsize='medium')

ax.figure.savefig('X_vs_Y-bar.png', dpi=400, bbox_inches='tight', pad_inches=0)
df.to_csv(f'{"X_vs_Y-out"}.csv', encoding='utf-8')
plt.show()
