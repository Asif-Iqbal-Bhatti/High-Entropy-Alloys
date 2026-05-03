#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert VASP results (via py4vasp vaspout.h5) to an Extended XYZ trajectory.

- Reads structure, forces, stress (and energy if available) in batches.
- Writes an .extxyz file compatible with ASE and other MD analysis tools.

Extended XYZ conventions used here:
  * Lattice="a11 a12 a13 a21 a22 a23 a31 a32 a33"  (Å)
  * Properties=species:S:1:pos:R:3:force:R:3
  * pbc="T T T"
  * Stress="s11 s12 s13 s21 s22 s23 s31 s32 s33" (units as provided by py4vasp)
  * Energy=<float> (if available from py4vasp)

Notes on units:
  - Positions are output in Å (fractional → cartesian via lattice).
  - Lattice is in Å.
  - Forces are written in eV/Å.
  - Stress comes from py4vasp; check your py4vasp/VASP version for units
    (commonly kBar). 

Author: Asif iqbal BHATTI/Asif_em2r
"""

def vaspouth5toextxyz():
    if os.path.exists(xyzpath):
        os.remove(xyzpath)    
    calc = Calculation.from_file("vaspout.h5")
    #print(dir(calc.energy[:]))
    #print(dir(calc))             
    #print(calc.structure[0].to_POSCAR())
    
    def to_extxyz(filename, calc, slice_):
        structure_data = calc.structure[slice_].read()
        forces = calc.force[slice_].read()["forces"]
        stress = calc.stress[slice_].read()["stress"]
    
        num_step, num_atom, _ = forces.shape
    
        species = [name.split("_")[0] for name in structure_data["names"]]
    
        with open(filename, "a") as xyz_file:
            for local_step in range(num_step):
                lattice = structure_data["lattice_vectors"][local_step]
                frac_positions = structure_data["positions"][local_step]
                positions = frac_positions @ lattice # MAT MUL in numpy
                step_forces = forces[local_step]
                step_stress = stress[local_step]
    
                xyz_file.write(f"{num_atom}\n")
    
                stress_flat = " ".join(f"{x:.6f}" for x in step_stress.flatten())
                lattice_flat = " ".join(f"{x:.6f}" for x in lattice.flatten())
                xyz_file.write(
                    f'Lattice="{lattice_flat}" '
                    f'Properties=species:S:1:pos:R:3:forces:R:3 '
                    f'pbc="T T T" '
                    f'Stress="{stress_flat}"\n'
                )
    
                for atom in range(num_atom):
                    sp = species[atom]
                    x, y, z = positions[atom]
                    fx, fy, fz = step_forces[atom]
                    xyz_file.write(f"{sp} {x:.6f} {y:.6f} {z:.6f} {fx:.6f} {fy:.6f} {fz:.6f}\n")
            
    # --- MAIN ---
    n_steps = calc.structure[:].number_steps()
    batch = 20000
    for start in range(0, n_steps, batch):
        stop = min(start + batch, n_steps)
        print(f"Writing steps {start}:{stop}")
        to_extxyz(EXTXYZ_NAME, calc, slice(start, stop))

vaspouth5toextxyz()
