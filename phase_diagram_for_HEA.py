#!/usr/bin/env python


import os, re, sys
import json
from pathlib import Path
from joblib import Parallel,delayed
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.phase_diagram import *
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.analysis.interface_reactions import InterfacialReactivity,GrandPotentialInterfacialReactivity
from pymatgen.io.vasp.outputs import Vasprun, Outcar 

initial_path = os.getcwd()
argo_path = 'alloy'
mpr = MPRester("YOUR_KEY")
compat 	= MaterialsProject2020Compatibility()
output = {}

def hea(path_main):
	path_main = Path(path_main)
	try:
		struct 	= Vasprun(os.path.join(path_main,'vasprun.xml'), parse_potcar_file = True)
		out_car	= os.path.abspath(os.path.join(path_main,'OUTCAR'))
	except: 
		print('---> File missing or corrupt: {}'.format(path_main))
		pass
		
	output['id'] = str(path_main).rsplit('/', 3)[1]
	g = struct.get_computed_entry(); 
	entry = compat.process_entries([g])[-1]
	entries = [entry]
	name = entry.name
	output['formula'] 										= entry.name
	output['natoms'] 											= round(entry.energy/entry.energy_per_atom)
	output['energy_per_atom'] 						= entry.energy_per_atom
	output['energy_per_atom_uncorrected'] = entry.uncorrected_energy/output['natoms']
	output['correction'] 									= output['energy_per_atom'] - output['energy_per_atom_uncorrected']
	
	line 			= re.findall('[A-Z][^A-Z]*', name.replace('(','').replace(')',''))
	searchset	= set(re.sub('\d',' ',' '.join(line)).split())
	temp 			= filter(lambda e: set(re.sub('\d',' ',str(e.composition).replace(' ','')).split())==searchset, entries)
	
	all_entries 		= mpr.get_entries_in_chemsys(set(searchset)) + list(temp)
	phase 					= PhaseDiagram(all_entries)
	print(phase)
	ehull 					= phase.get_e_above_hull(entry) 
	output['Ehull'] = round(ehull, 5); print(ehull) # eV/atom
	return output
	
if __name__ == '__main__':
	abs_path = os.path.abspath(os.path.join(initial_path, argo_path))
	argo_dir = next(os.walk(abs_path))[1]
	def runthis(p):
		arg_path = Path(os.path.join(argo_path, p))
		print(arg_path)
		try:
			out = hea(str(arg_path))
		except:
			out = dict({'fail_': 'fail_' + str(p)})
		return out
	out = Parallel(n_jobs = -1, backend="threads", batch_size='auto')(delayed(runthis)(p) for p in argo_dir)
	json.dump(out, open(argo_path + '.json','w'), indent = 1, sort_keys=False, ensure_ascii=True)

