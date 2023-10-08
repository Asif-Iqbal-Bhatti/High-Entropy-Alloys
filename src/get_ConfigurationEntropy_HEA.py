#!/usr/bin/env python3

'''
#============================================================================
# AUTHOR:: Asif IQBAL
# GitHub:: @AIB_EM
# USAGE :: TO compute the partial or compositional configurational entropic
#          term. The definition is ∆Sconf = -kb * SUMi (ci*ln(ci)) [eV/K]
#============================================================================
'''

from ase.io import read
import pprint
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ase.units import kB,mol,kJ

print(f"VALUE OF kB: {kB:15.8e}") #kB = 8.617333262E-5 #eV/K  # Boltzmann constant

def get_ENTROPY_from_CIF():
    cifFile = '../LPSC_cif_files/2017/Li6PS5Cl_structure_from_2017_paper.cif'
    structureList = read(cifFile, format = 'cif', index = ':', store_tags = True, reader = 'ase')[-1] #, fractional_occupancies = True
    
    cifTags = structureList.info.copy()
    
    siteSymbols    = np.array(cifTags['_atom_site_type_symbol'])
    siteOccupancy = np.array(cifTags['_atom_site_occupancy'])
    siteMulti       = np.array(cifTags['_atom_site_symmetry_multiplicity'])
    siteWyckoff       = np.array(cifTags['_atom_site_wyckoff_symbol'])
    
    result_dict = {}
    for i, symbol in enumerate(siteSymbols):
        occupancy = siteOccupancy[i]
        multi = siteMulti[i]
        wyckoff = siteWyckoff[i]
        
        key = f'{symbol}_{multi}{wyckoff}'
        value = [occupancy, multi, wyckoff, f"{multi}{wyckoff}"]
      
        result_dict[key] = value
        
    pp = pprint.PrettyPrinter(indent=4)    
    pp.pprint(result_dict)

    # Create a dictionary to filter keys based on Wycoff Multiplicity
    third_values_dict = {}
    
    for key, value in result_dict.items():
        third_value = value[2]
        
        if third_value not in third_values_dict:
            third_values_dict[third_value] = [key]
        else:
            third_values_dict[third_value].append(key)
    
    common_Wyckoff = [keys for keys in third_values_dict.values() if len(keys) > 1]
    print(f"{'Same Wycoff positions':15s}: {common_Wyckoff}")
    
    # sublattice configuration entropy approximation
    # Sum over sublattice index with conentration
    ee = 0.0
    for k, v in result_dict.items():
        if k in np.array(common_Wyckoff)[:,0]:
            x = v[0]
            ee += v[1] * ( x * np.log(x) + (1-x) * np.log(1-x) )
    oo = -kB * ee
    
    r1 = []
    for t in np.linspace(0, 1000, num=101):
        ee_per_atom = (oo * t) / 52
        r1.append((t, ee_per_atom))
    
    df1 = pd.DataFrame(r1, columns=['T[K]', 'TS[eV/at]'])
    df1.to_csv("LPSCCONFENTROPY.csv", index=False)
    df1.plot(kind='scatter', x='T[K]', y='TS[eV/at]', color='red')
    plt.show()

get_ENTROPY_from_CIF()

    
