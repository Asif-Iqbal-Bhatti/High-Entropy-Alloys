#!/usr/bin/env python3

"""
Madelung Potential Calculator for VASP structures.

This script parses a VASP POSCAR/CONTCAR file, extracts atomic positions 
and charges (from the 4th column), and computes the Madelung potential 
at each atomic site using a direct-sum supercell approach.

Author: Asif Iqbal
"""

import numpy as np
import os
import sys

# Physical Constants (CODATA 2018)
EPS0 = 8.8541878128e-12    
E_CHG = 1.602176634e-19    
A2M = 1e-10                
KE = 1 / (4 * np.pi * EPS0) 

def read_vasp(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Structure file '{filename}' not found.")
    with open(filename, 'r') as f:
        lines = f.readlines()
    try:
        scale = float(lines[1].strip())
        lattice = np.array([l.split() for l in lines[2:5]], dtype=float) * scale
        counts = list(map(int, lines[6].split()))
        total_atoms = sum(counts)
        offset = 8 if "selective" in lines[7].lower() else 7
        coord_type = lines[offset].strip().lower()
        coords, charges = [], []
        for i in range(total_atoms):
            parts = lines[offset + 1 + i].split()
            coords.append(parts[:3])
            charges.append(float(parts[3]) if len(parts) >= 4 else 0.0)
        coords = np.array(coords, dtype=float)
        charges = np.array(charges, dtype=float)
        if coord_type.startswith('d'):
            coords = coords @ lattice
        return lattice, coords, charges
    except Exception as e:
        raise ValueError(f"Failed to parse VASP file: {e}")

def calculate_madelung(lattice, coords, charges, shell_limit=3):
    n_atoms = len(coords)
    potentials = np.zeros(n_atoms)
    net_charge = np.sum(charges)
    if abs(net_charge) > 1e-3:
        print(f"Warning: System not neutral (Net: {net_charge:.3f}).")
    r_ranges = range(-shell_limit, shell_limit + 1)
    offsets = np.array([i*lattice[0] + j*lattice[1] + k*lattice[2] 
               for i in r_ranges for j in r_ranges for k in r_ranges])
    for i in range(n_atoms):
        phi_i = 0.0
        pos_i = coords[i]
        for image_offset in offsets:
            pos_j_images = coords + image_offset
            dists = np.linalg.norm(pos_j_images - pos_i, axis=1)
            for j in range(n_atoms):
                if dists[j] < 1e-5: continue
                phi_i += charges[j] / dists[j]
        potentials[i] = phi_i
    return potentials * ((E_CHG / A2M) * KE)

def main():
    input_file = 'CONTCAR' if os.path.exists('CONTCAR') else 'NaCl.poscar'
    try:
        lattice, coords, charges = read_vasp(input_file)
        print(f"Reading {input_file}: {len(coords)} atoms found.")
        potentials = calculate_madelung(lattice, coords, charges)
        total_energy_ev = 0.5 * np.sum(charges * potentials)
        print("\n" + "="*50)
        print(f"{'Index':<6} {'Charge (e)':<12} {'Potential (V)':<15}")
        print("-"*50)
        for i, (q, v) in enumerate(zip(charges, potentials)):
            print(f"{i+1:<6} {q:<12.3f} {v:<15.4f}")
        print("="*50)
        print(f"Total Madelung Energy: {total_energy_ev:.6f} eV")
        print("="*50)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
