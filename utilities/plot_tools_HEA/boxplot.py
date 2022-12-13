#!/usr/bin/env python3

'''
#=========================================================
# AUTHOR:: ASIF -> AIB_EM
# USAGE :: PLOT BOX PLOT  
#=========================================================
'''

import pandas as pd, numpy as np, seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

maindata = r'Sheet1.csv'
sns.set_style('white')
sns.set_theme(style="ticks")
sns.set_context("paper", font_scale = 1.1)
#print(plt.style.available)
#plt.style.use('bmh') 
plt.rcParams.update({
	"text.usetex": True,
	"font.family": "Times",
	"font.sans-serif": ['Helvetica', 'DejaVu Sans', 'Verdana'],
})

gg = []
df = pd.read_csv(maindata)
hdd = df.columns.values.tolist()
for i in (hdd):	
	gg.append( list(df[i].dropna()) )

def statplot_cols():
	xnames=[
	'(6,0,0)','(6,1,0)','(6,1,1)','(8,0,0)',
	'(6,0,0)','(6,1,0)','(6,1,1)','(8,0,0)',
	'(6,0,0)','(6,1,0)','(6,1,1)','(8,0,0)',
	'(6,0,0)','(6,1,0)','(6,1,1)','(8,0,0)','(8,1,0)','(8,1,0)',
	'(6,0,0)','(6,1,0)','(6,1,1)','(8,0,0)',
	'(6,0,0)','(6,1,0)','(6,1,1)','(8,0,0)'
	]
	
	pp = sns.boxplot(data=df.values, showfliers=False, showmeans=True, 
							meanprops={"marker":"o",
							"markerfacecolor":"white",
							"markeredgecolor":"black",
							"markersize":"0.5"},
							dodge=True, \
	color='darkgrey', whis=np.inf, width=0.4, linewidth=0.2)
	
	means1 = [np.mean(xi) for xi in gg]
	plt.plot(hdd, means1, '--', lw=0.5, color='blue')
	
	max1 = [np.max(xi) for xi in gg]
	plt.plot(hdd, max1, 'o', lw=0.5, color='black', markersize=4)
	
	min2 = [np.min(xi) for xi in gg]
	plt.plot(hdd, min2, 'o-', mfc='none', lw=0.5, color='green', markersize=4)
	
	#==============================================================
	# x= 0.0
	left, bottom, width, height = (-0.2, 74, 3.5, 50)
	rect=mpatches.Rectangle((left,bottom),width,height,alpha=0.1, facecolor="grey")
	plt.gca().add_patch(rect)
	plt.text(0, 120,'${x=0.0}$',fontsize=12, color="black", weight="bold")
	
	# x= 0.11
	left, bottom, width, height = (3.7, 74, 3.6, 50)
	rect=mpatches.Rectangle((left,bottom),width,height,alpha=0.1, facecolor="grey")
	plt.gca().add_patch(rect)
	plt.text(3.8, 120,'${x=0.11}$',fontsize=12, color="black", weight="bold")
	
	# x= 0.20
	left, bottom, width, height = (7.7, 74, 3.6, 50)
	rect=mpatches.Rectangle((left,bottom),width,height,alpha=0.1, facecolor="grey")
	plt.gca().add_patch(rect)
	plt.text(7.9, 120,'${x=0.20}$',fontsize=12, color="black", weight="bold")
	
	# x= 0.33
	left, bottom, width, height = (11.7, 74, 5.6, 50)
	rect=mpatches.Rectangle((left,bottom),width,height,alpha=0.1, facecolor="grey")
	plt.gca().add_patch(rect)
	plt.text(13, 120,'${x=0.33}$',fontsize=12, color="black", weight="bold")
	
	# x= 0.40
	left, bottom, width, height = (17.7, 74, 3.6, 50)
	rect=mpatches.Rectangle((left,bottom),width,height,alpha=0.1, facecolor="grey")
	plt.gca().add_patch(rect)
	plt.text(18, 120,'${x=0.40}$',fontsize=12, color="black", weight="bold")
	
	# x= 0.50
	left, bottom, width, height = (21.6, 74, 3.7, 50)
	rect=mpatches.Rectangle((left,bottom),width,height,alpha=0.1, facecolor="grey")
	plt.gca().add_patch(rect)
	plt.text(22, 120,'${x=0.50}$',fontsize=12, color="black", weight="bold")
	
	pp.set_xticklabels(xnames, rotation=90)
	pp.set_xlabel('(p,t,q)',fontsize=12)
	pp.set_ylabel('H[meV/atom]',fontsize=12)
	plt.savefig('test.png', dpi=600,bbox_inches='tight')

#==================================================
#    DEF 2
#==================================================
	
def boxplot_2():
	ff = {
	'0.0': {
	'(6,0,0)0':gg[0],
	'(6,1,0)0':gg[1],
	'(6,1,1)0':gg[2],
	'(8,0,0)0':gg[3]},
	'0.11':{
	'(6,0,0)1':gg[4],
	'(6,1,0)1':gg[5],
	'(6,1,1)1':gg[6],
	'(8,0,0)1':gg[7]},
	'0.20': {
	'(6,0,0)2':gg[8],
	'(6,1,0)2':gg[9],
	'(6,1,1)2':gg[10],
	'(8,0,0)2':gg[11]},
	'0.33': {
	'(6,0,0)3':gg[12],
	'(6,1,0)3':gg[13],
	'(6,1,1)3':gg[14],
	'(8,0,0)3':gg[15],
	'(8,1,0)3':gg[16],
	'(8,1,1)3':gg[17]},
	'0.40': {
	'(6,0,0)4':gg[18],
	'(6,1,0)4':gg[19],
	'(6,1,1)4':gg[20],
	'(8,0,0)4':gg[21]},
	'0.50': {
	'(6,0,0)5':gg[22],
	'(6,1,0)5':gg[23],
	'(6,1,1)5':gg[24],
	'(8,0,0)5':gg[25]},
	}
	
	plt.clf()
	df = pd.DataFrame()
	appended_data=[]
	for comp, ptq_groups in ff.items():
		for ptq, values in ptq_groups.items():
			data = pd.DataFrame({'H [meV/atom]': values, 'ptq': ptq, 'comp': comp}); 
			appended_data.append(data)
	df = pd.concat(appended_data)
	g=sns.boxplot(data=df, y='ptq', x='H [meV/atom]', palette='spring', showfliers=False)
	
	plt.savefig('test1.png', dpi=300,bbox_inches='tight')



if __name__ == '__main__':
	statplot_cols()
