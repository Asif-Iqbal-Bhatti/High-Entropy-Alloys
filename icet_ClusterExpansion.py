#!/usr/bin/env python3
	
import os
from ase.db import connect
from ase.io import read
from tqdm import tqdm
import numpy as np 

'''
cwd = os.path.abspath(os.path.dirname(__file__))
outdbpathF = os.path.join(cwd, 'reference_data.db')

#============= CREATE DATABASE
inp_Argo = 'VASPRun_MAPS'
initial_path = os.getcwd()
abs_path = os.path.abspath(os.path.join(initial_path, inp_Argo))
argo_Dir = next(os.walk(abs_path))[1]
print(f'REF DATA:: {str(len(argo_Dir))}')

d=0
with connect(outdbpathF) as conDB:
	for p in tqdm(argo_Dir):
		arg_path  = os.path.join(inp_Argo, p)
		vasp_xml	= os.path.join(arg_path,'vasprun.xml')
		read_xmlF = read(vasp_xml, index='-1') # FINAL/ALL GEOMTERY (USE "-1" or ":")	
		
		g = read_xmlF.get_cell()
		lp = [np.linalg.norm(g[0]),np.linalg.norm(g[1]),np.linalg.norm(g[2])]
		
		conDB.write(read_xmlF, \
		tag = read_xmlF.get_chemical_formula() + '_' + str(d), \
		lattice_parameter = np.mean(lp), \
		tot_energy = read_xmlF.get_total_energy()
		)			
		
		d+=1
'''

from icet import ClusterSpace, StructureContainer, ClusterExpansion
from trainstation import CrossValidationEstimator

db = connect('reference_data.db')
primitive_structure = db.get(id=10).toatoms()  # primitive structure
cs = ClusterSpace(structure=primitive_structure,
                  cutoffs=[13, 6.5, 6],
                  chemical_symbols=['A', 'B', 'C', 'D', 'E'])


sc = StructureContainer(cs)
for row in db.select():
	sc.add_structure(structure=row.toatoms(),
                    user_tag=row.tag,
                    properties={'tot_energy': row.tot_energy})
	print(sc)
print(sc)
