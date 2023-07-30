#!/usr/bin/env python3

'''
#==================================================
# USAGE  :: Regression analysis 
# AUTHOR :: ASIF IQBAL
# DATED  :: 12/7/2023
#==================================================
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


KB = 8.617333262E-5 #eV K-1 

outFILE = 'Tix0.0_SRO_vs_TE_'
inpFILE = 'VASPRun_x_0.00_800_1st_2nd_3rd.csv'
df = pd.read_csv(inpFILE)
orig = df.iloc[:, -3:]
df_s = orig.round(decimals=5)
unique_df = df_s.drop_duplicates(subset=['Total_SRO', 'TE/at'])
orig_sor = unique_df.sort_values(by=['TE/at'], ascending=True)
x_data = (orig_sor["TE/at"] - orig_sor["TE/at"].max()) * 1000
y_data = orig_sor["Total_SRO"]

def plot_SRO_vsd_Energy():
    fig, ax = plt.subplots(figsize=(8,6))
    print(len(x_data), len(y_data))
    max_index = orig_sor['Total_SRO'].idxmax()
    min_index = orig_sor['Total_SRO'].idxmin()
    max_point = orig_sor.loc[max_index]
    min_point = orig_sor.loc[min_index]
    
    plt.scatter(max_point['TE/at'], max_point['Total_SRO'], color='red', label='Max Point')
    plt.scatter(min_point['TE/at'], min_point['Total_SRO'], color='blue', label='Min Point')
    # ------
    coefficients = np.polyfit(x_data, y_data, deg=1)
    fit_line = np.polyval(coefficients, x_data)
    plt.plot(x_data, fit_line, color='red', label='Fitted Line')
    # ------
    
    plt.scatter(x_data, y_data, c='k', label='Data', facecolors='none', edgecolors='magenta', linewidths=1,)
    plt.xlabel('Relative TE [$meV/atom$]')
    plt.ylabel('SRO')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator(10))
    plt.ylim(-0.1,2.0)
    plt.xlim(-15, 1)
    plt.legend()
    plt.grid(False)
    plt.title('SRO vs Te/atom for $Ti_{x} = 0.0$') 
    plt.savefig(outFILE+"E.pdf", dpi=300,bbox_inches="tight")
    plt.show()


def umap_FIT():
    fig, ax = plt.subplots(figsize=(8,6))
    # UMAP FITTING
    import umap
    data = np.column_stack((x_data, y_data))
    umap_embedding = umap.UMAP(n_components=2, random_state=42).fit_transform(data)
    plt.scatter(umap_embedding[:, 0], umap_embedding[:, 1])
    plt.xlabel('UMAP Dimension 1')
    plt.ylabel('UMAP Dimension 2')
    plt.title('UMAP Projection of Scatter Plot Data')
    plt.grid(True)
    plt.savefig(outFILE+"UMAP.pdf", dpi=300,bbox_inches="tight")
    plt.clf()


def gpr_FIT():
    fig, ax = plt.subplots(figsize=(8,6))
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
    
    X = x_data.values.reshape(-1, 1)  # GPR requires 2D input
    y = y_data.values
    
    # Define the Gaussian Process kernel (RBF kernel with automatic hyperparameter optimization)
    kernel = C(1.0, (1e-3, 1e3)) * RBF(1.0, (1e-2, 1e2))
    
    # Fit Gaussian Process Regression
    gpr = GaussianProcessRegressor(
                    kernel=kernel, 
                    n_restarts_optimizer=20, 
                    random_state=42,
                    optimizer = 'fmin_l_bfgs_b').fit(X, y)
    
    # Generate predictions and uncertainty estimates
    x_pred = np.linspace(min(X), max(X), 100).reshape(-1, 1)
    y_pred, y_std = gpr.predict(x_pred, return_std=True)
    
    # Plot the scatter plot and the fitted curve with uncertainty bounds
    plt.scatter(X, y, label='Data', facecolors='none', edgecolors='k', linewidths=1,)
    plt.plot(x_pred, y_pred, color='red', label='GPR Fitted Curve')
    plt.fill_between(x_pred.flatten(), y_pred - 1.96 * y_std, y_pred + 1.96 * y_std, color='gray', alpha=0.2, label='Uncertainty Bounds')
    plt.xlabel('TE/at')
    plt.ylabel('Total_SRO')
    plt.legend()
    plt.grid(True)
    plt.savefig(outFILE+"GPR.pdf", dpi=300,bbox_inches="tight")
    plt.show()
    plt.clf()


'''
# METROPLOIS MONTE CARLO ALGORITHM
# 
#
'''

def metropolis_monte_carlo():
    fig, ax = plt.subplots(figsize=(8,6))
    df = pd.read_csv(inpFILE)
    orig = df.iloc[:, -3:]
    df_s = orig.round(decimals=5)
    unique_df = df_s.drop_duplicates(subset=['Total_SRO', 'TE/at'])
    orig_sor = unique_df.sort_values(by=['TE/at'], ascending=True)
    
    energies = orig_sor['TE/at']
    old_energy = np.random.choice(energies)
    #old_energy = -10.735077
    print("REF ENE:", old_energy)
    accepted_energies = []

    temp = 0.1
    for _ in range(int(1E9)):
        new_energy = np.random.choice(energies)
        x = np.exp(-(new_energy - old_energy) / (KB*temp))
        
        if new_energy <= old_energy:
            #ind1 = orig[orig.isin([new_energy]).any(axis=1)]
            ind1 = orig_sor[orig_sor['TE/at'] == new_energy]
            #print(ind1)
            accepted_energies.append({'TE/at': ind1['TE/at'].values[0],'Total_SRO': ind1['Total_SRO'].values[0],})
            old_energy = new_energy
        elif x >= np.random.rand():
            ind1 = orig_sor[orig_sor['TE/at'] == new_energy]
            accepted_energies.append({'TE/at': ind1['TE/at'].values[0],'Total_SRO': ind1['Total_SRO'].values[0],})
            old_energy = new_energy

            
    df1 = pd.DataFrame(accepted_energies)
    df1 = df1.drop_duplicates()
    df1.to_csv('TEST.csv', index=False)
    print(df1)    
    x_data = df1['TE/at']
    y_data = df1['Total_SRO']
    coefficients = np.polyfit(x_data, y_data, deg=1)
    fit_line = np.polyval(coefficients, x_data)
    plt.plot(x_data, fit_line, color='red', label='Fitted Line')
    plt.scatter(x_data, y_data, c='k', label='Data', facecolors='none', edgecolors='green', linewidths=1,)
    plt.xlabel('Relative TE [$eV/atom$]')
    plt.ylabel('SRO')
    plt.legend()
    plt.grid(False)
    plt.title('SRO vs Te/atom for $Ti_{x} = 0.0$') 
    plt.savefig(outFILE+"MC.pdf", dpi=300,bbox_inches="tight")
    plt.show()   


plot_SRO_vsd_Energy()
#umap_FIT()
#gpr_FIT()
#metropolis_monte_carlo()



