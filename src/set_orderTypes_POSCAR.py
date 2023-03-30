#!/usr/bin/env python3


from ase.io import read, write

atoms = read('POSCAR')
symbols = atoms.get_chemical_symbols()
new_order = ['Cl', 'P', 'S', 'Li'] # new order of atom types

new_indices = []

for symbol in new_order:
    new_indices += [i for i, s in enumerate(symbols) if s == symbol]

print(new_indices)
atoms = atoms[new_indices]

write('POSCAR_reordered', atoms, format='vasp', direct=True)
