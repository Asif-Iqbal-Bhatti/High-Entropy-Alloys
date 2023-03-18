#!/usr/bin/env python3

import matplotlib.pyplot as plt
import seaborn as sns

   
def plot_POSCAR_using_pmg():
    from pymatgen.core.structure import Molecule, Structure
    from pymatviz import plot_structure_2d
    from pymatgen.io.vasp import Vasprun, Outcar, Poscar
    fig = plt.figure()
    fig, ax = plt.subplots(1, 2, figsize=(8, 6))  
    struct = Structure.from_file('POSCAR_grid')
    tit = Poscar.from_file('POSCAR_grid', check_for_POTCAR=False).comment
    ll = struct.lattice.lengths
    a, b = ll[0], ll[1]
    
    plot_structure_2d(struct, ax[0], rotation='0x,0y,0z',
                    atomic_radii={'H':0.05, 'Fe':0.2}, 
                    colors={'H':'gray','Fe':'b'},
                    show_unit_cell=True, site_labels=False,
                    label_kwargs={"fontsize": 4,})
                    
    ax[0].set_title(f"Grid sampling of Fe-bcc {tit}")
    ax[0].axhline(y=a/2, color='black', linestyle='--', linewidth=0.2)
    ax[0].axvline(x=b/2, color='black', linestyle='--', linewidth=0.2)
    ax[0].set_xlabel("x-axis [$\mathrm{\AA}$]")
    ax[0].set_ylabel("y-axis [$\mathrm{\AA}$]")
    ax[0].set_xlim(0, a)
    ax[0].set_ylim(0, b)

    plot_structure_2d(struct, ax[1], rotation='0x,90y,0z',
                    atomic_radii={'H':0.05, 'Fe':0.2}, 
                    colors={'H':'gray','Fe':'g'},
                    show_unit_cell=True, site_labels=False,
                    label_kwargs={"fontsize": 8})

    plt.savefig(f"h_grid_sample.pdf", bbox_inches='tight', dpi=300)

plot_POSCAR_using_pmg()

def plot_POSCAR_using_ASE():
    from ase.io import read
    
    structure = read('POSCAR_grid')
    cell_size = structure.get_cell_lengths_and_angles()[:3]
    coordinates = structure.get_positions()
    atom_types = structure.get_chemical_symbols()
    atom_props = {
        'H': {'color': 'red', 'size': 1.5},
        'Fe': {'color': 'blue', 'size': 2.0},
    }
    
    fig, ax = plt.subplots()
    for i, atom_type in enumerate(atom_types):
        color = atom_props[atom_type]['color']
        x, y, z = coordinates[i]
        ax.scatter(x, y, s=30.5, c=color)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Atomic Coordinates')
    
    ax.set_aspect('equal')
    ax.set_xlim([0, cell_size[0]])
    ax.set_ylim([0, cell_size[1]])
    plt.show()
