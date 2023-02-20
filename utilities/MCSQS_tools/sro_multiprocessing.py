#!/usr/bin/env python
#
# USING pyscal to compute SRO parameters using fontain et al.
# https://pyscal.org/en/latest/_modules/pyscal/core.html#System.calculate_pmsro
#
import sys, os
import numpy as np
import pyscal.core as pc
from ase.io.vasp import write_vasp, read_vasp
from tqdm import tqdm
import concurrent.futures

top = pc.System()
cwd = os.path.abspath(os.path.dirname(__file__)) 
directories = [os.path.join(cwd, k) for k in os.listdir('.') if os.path.isdir(k)]

def calculate_sro(pair):
    t, j = pair
    sro = top.calculate_pmsro(reference_type=t+1, compare_type=j+1, average=True, shells=2, delta=True)
    pair_name = f"{atom_type[t]}-{atom_type[j]}"
    sro_1st_shell = round(sro[0], 4)
    sro_2st_shell = round(sro[1], 4)
    sro_1st_2nd_shell = round(sro[0]+sro[1], 4)
    return [pair_name, sro_1st_shell, sro_2st_shell, sro_1st_2nd_shell]

with open("SRO_column.dat", "w") as f:
    for k in tqdm(next(os.walk('.'))[1]):
        d = os.path.join(cwd, k)

        top.read_inputfile(os.path.join(d, 'POSCAR_bestsqs'), format='poscar')
        top.find_neighbors(method='cutoff', cutoff=1.2)
        
        initial = read_vasp(os.path.join(d, 'POSCAR_bestsqs'))
        at_type = initial.get_chemical_symbols()
        atom_type = sorted(list(set(at_type)))
        
        pairs_data = []
        header = "{:6s}  {:10.10s} {:10.10s} {:10.10s}\n".format("Pairs", "1st Shell", "2nd Shell", "1st+2nd")
        f.write("{}\n".format(k))
        f.write(header)
        
        pairs = [(t, j) for t in range(5) for j in range(t, 5)]
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = list(executor.map(calculate_sro, pairs))
        pairs_data.extend(results)
        
        for data in pairs_data:
            format_str = "{} {:9.5f} {:9.5f} {:9.5f}\n"
            f.write(format_str.format(data[0], data[1], data[2], data[3]))
                    
        a = sum([data[1] for data in pairs_data])
        b = sum([data[2] for data in pairs_data])
        c = sum([data[3] for data in pairs_data])
    
        f.write("{}\n".format("-"*35))
        f.write ("{:5s} {:9.5f} {:9.5f} {:9.5f}\n".format('T_SRO', a, b, a+b))
        f.write("{}\n".format("-"*50))
