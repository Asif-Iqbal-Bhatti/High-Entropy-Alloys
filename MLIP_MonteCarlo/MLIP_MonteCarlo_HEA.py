#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Monte Carlo Simulation for High Entropy Structure Optimization
(using ASE + MLIP instead of VASP)

- Generates substituted supercell configurations (high-entropy)
- Evaluates energies using an ASE MLIP calculator (GRACE/TensorPotential via TPCalculator)
- Uses Metropolis Monte Carlo acceptance at temperature T
- Avoids duplicate configurations using a structure fingerprint
- Logs accepted/rejected energies and writes POSCAR snapshots

Original Author: Noah Oyeniran
Version: 1:0:0
License: MIT

FURTHER MODIFICATION by Author :: Asif_EM2R
"""

import os
import random
import hashlib
import numpy as np

import matplotlib
matplotlib.use("Agg")  # must be BEFORE importing pyplot
import matplotlib.pyplot as plt


from pymatgen.core import Structure, Element
from pymatgen.io.vasp import Poscar
from pymatgen.io.ase import AseAtomsAdaptor

from ase.constraints import FixAtoms
from ase.optimize import FIRE2

# === GRACE / TensorPotential ASE calculator ===
from tensorpotential.calculator.foundation_models import grace_fm

# ---------------------------
# User configuration
# ---------------------------

MC_TEMPERATURE = 300.0  # Kelvin
BOLTZMANN_CONSTANT = 8.617333262145e-5  # eV/K
MAX_ITERATIONS = 5000

# If True, do a short relaxation (positions only) before energy evaluation
RELAX_BEFORE_ENERGY = False
RELAX_FMAX = 0.05   # eV/Å
RELAX_STEPS = 200

# Neighbor cutoff for fingerprint (Å)
FINGERPRINT_R = 3.0

# RNG seed (optional)
RANDOM_SEED = 42

# ---------------------------
# Helpers
# ---------------------------

def reorder_species_concisely(structure: Structure) -> Structure:
    """Reorders species in a structure for consistent POSCAR output."""
    composition = structure.composition
    sorted_elements = sorted(composition.as_dict().keys())
    reordered_sites = []
    for el in sorted_elements:
        el_sites = [site for site in structure if site.species_string == el]
        reordered_sites.extend(el_sites)
    return Structure.from_sites(reordered_sites)


def _compute_integer_counts(total_sites: int, fractions_dict: dict) -> dict:
    """
    Convert mole fractions to integer substitution counts that sum <= total_sites.
    Uses floor + remainder distribution by largest fractional parts.
    """
    # raw counts
    raw = {el: fractions_dict[el] * total_sites for el in fractions_dict}
    base = {el: int(np.floor(raw[el])) for el in raw}

    used = sum(base.values())
    remaining = max(0, min(total_sites - used, total_sites))

    # distribute remaining by fractional parts (largest first)
    frac_parts = sorted(((el, raw[el] - base[el]) for el in raw),
                        key=lambda x: x[1], reverse=True)

    i = 0
    while remaining > 0 and i < len(frac_parts):
        el, _ = frac_parts[i]
        base[el] += 1
        remaining -= 1
        i += 1

    # If fractions sum < 1.0, leftover sites stay as original species (not substituted).
    # If fractions sum > 1.0, the above caps by total_sites, so it won't exceed total_sites.
    return base


def create_supercell_and_substitute(
    poscar_file: str,
    scaling_matrix,
    substitutions: dict,
    x_elements: list,
    output_dir: str
    ) -> Structure:
    """
    Creates a supercell and performs element substitutions.

    substitutions example:
      {"Mo": {"Cr": 0.33, "Ti": 0.33}}

    x_elements:
      species that should not be substituted (fixed sublattices, terminations)
    """
    structure = Structure.from_file(poscar_file)
    structure.make_supercell(scaling_matrix)

    x_elements = {Element(el) for el in x_elements}

    for element_to_replace, substituent_data in substitutions.items():
        replaceable_sites = [
            i for i, site in enumerate(structure)
            if site.species_string == element_to_replace
            and Element(site.species_string) not in x_elements
        ]

        total_sites = len(replaceable_sites)
        if total_sites == 0:
            continue

        random.shuffle(replaceable_sites)

        # Determine integer substitution counts
        counts = _compute_integer_counts(total_sites, substituent_data)

        start = 0
        for sub_el, nsub in counts.items():
            if nsub <= 0:
                continue
            idxs = replaceable_sites[start:start + nsub]
            for site_index in idxs:
                structure[site_index] = Element(sub_el)
            start += nsub

    structure = reorder_species_concisely(structure)

    # write current POSCAR (for inspection)
    os.makedirs(output_dir, exist_ok=True)
    Poscar(structure).write_file(os.path.join(output_dir, "POSCAR"))

    return structure


def generate_structure_fingerprint(structure: Structure, r: float = 3.0) -> str:
    """
    Generates a heuristic fingerprint:
      species + coordination number + neighbor species counts within radius r
    """
    site_data = []
    for site in structure:
        neighbors = structure.get_neighbors(site, r=r)
        coord = len(neighbors)
        species_count = {}
        for nbr in neighbors:
            s = nbr.species_string
            species_count[s] = species_count.get(s, 0) + 1
        site_data.append((site.species_string, coord, tuple(sorted(species_count.items()))))

    return hashlib.md5(str(site_data).encode("utf-8")).hexdigest()


def get_mlip_calculator() -> grace_fm:
    model_path = 'GRACE-2L-OAM'
    calc = grace_fm(
        model_path, 
        pad_atoms_number=2, 
        pad_neighbors_fraction=0.05, 
        min_dist=0.5, 
        device='cuda'  # 'cuda' for GPU, 'cpu' for CPU
    )
    return calc


def energy_from_mlip(
    structure: Structure,
    calc,
    relax: bool = False,
    x_elements: list = None,
    fmax: float = 0.05,
    steps: int = 200
    ) -> tuple[float, Structure]:
    """
    Compute energy (eV) from MLIP, optionally relaxing positions with ASE.
    ASE workflow: attach calculator -> get_potential_energy. [1](https://ase-lib.org/ase/calculators/calculators.html)
    """
    adaptor = AseAtomsAdaptor()
    atoms = adaptor.get_atoms(structure)

    atoms.calc = calc

    if relax:
        if x_elements:
            fixed = [i for i, sym in enumerate(atoms.get_chemical_symbols()) if sym in set(x_elements)]
            if fixed:
                atoms.set_constraint(FixAtoms(indices=fixed))

        opt = FIRE(atoms, logfile=None)
        opt.run(fmax=fmax, steps=steps)

    e = atoms.get_potential_energy()  # eV
    new_structure = adaptor.get_structure(atoms)
    return float(e), new_structure


# ---------------------------
# Main Monte Carlo
# ---------------------------

def monte_carlo_mlip(
    poscar_file: str,
    scaling_matrix,
    substitutions: dict,
    x_elements: list,
    output_dir: str
    ):
    if RANDOM_SEED is not None:
        random.seed(RANDOM_SEED)
        np.random.seed(RANDOM_SEED)

    os.makedirs(output_dir, exist_ok=True)
    accepted_dir = os.path.join(output_dir, "accepted")
    rejected_dir = os.path.join(output_dir, "rejected")
    os.makedirs(accepted_dir, exist_ok=True)
    os.makedirs(rejected_dir, exist_ok=True)

    accepted_log = os.path.join(output_dir, "accepted_energies.txt")
    rejected_log = os.path.join(output_dir, "rejected_energies.txt")

    with open(accepted_log, "w") as f:
        f.write("Step\tEnergy(eV)\tDeltaE(eV)\tAcceptedProb\n")
    with open(rejected_log, "w") as f:
        f.write("Step\tEnergy(eV)\tDeltaE(eV)\tAcceptedProb\n")

    # Build calculator once
    calc = get_mlip_calculator()  
    
    # Initial structure
    current_structure = create_supercell_and_substitute(
        poscar_file, scaling_matrix, substitutions, x_elements, output_dir
    )

    unique = set()

    # Ensure initial structure is unique
    fp0 = generate_structure_fingerprint(current_structure, r=FINGERPRINT_R)
    unique.add(fp0)

    # Initial energy
    current_energy, current_structure = energy_from_mlip(
        current_structure,
        calc=calc,
        relax=RELAX_BEFORE_ENERGY,
        x_elements=x_elements,
        fmax=RELAX_FMAX,
        steps=RELAX_STEPS
    )

    best_energy = current_energy
    best_structure = current_structure

    energies_accepted = [current_energy]

    # Save initial
    Poscar(current_structure).write_file(os.path.join(accepted_dir, "POSCAR_step_0000.vasp"))
    Poscar(best_structure).write_file(os.path.join(output_dir, "BEST_POSCAR.vasp"))

    print(f"Initial energy: {current_energy:.6f} eV")

    kT = BOLTZMANN_CONSTANT * MC_TEMPERATURE

    for step in range(1, MAX_ITERATIONS + 1):
        print(f"\nStep {step}/{MAX_ITERATIONS}")

        # Propose a new candidate by regenerating random substitutions
        candidate_structure = create_supercell_and_substitute(
            poscar_file, scaling_matrix, substitutions, x_elements, output_dir
        )

        # Avoid duplicates (try a few times)
        tries = 0
        fp = generate_structure_fingerprint(candidate_structure, r=FINGERPRINT_R)
        while fp in unique and tries < 20:
            tries += 1
            candidate_structure = create_supercell_and_substitute(
                poscar_file, scaling_matrix, substitutions, x_elements, output_dir
            )
            fp = generate_structure_fingerprint(candidate_structure, r=FINGERPRINT_R)

        if fp in unique:
            print("Could not find a new unique candidate after 20 retries; skipping.")
            continue

        # Evaluate energy with MLIP
        try:
            cand_energy, cand_structure_relaxed = energy_from_mlip(
                candidate_structure,
                calc=calc,
                relax=RELAX_BEFORE_ENERGY,
                x_elements=x_elements,
                fmax=RELAX_FMAX,
                steps=RELAX_STEPS
            )
        except Exception as e:
            print(f"MLIP evaluation failed: {e}")
            continue

        dE = cand_energy - current_energy
        accept_prob = 1.0 if dE <= 0 else float(np.exp(-dE / kT))

        # Metropolis decision
        u = random.random()
        accepted = (u <= accept_prob)

        # Record uniqueness only after evaluation (so we don't block)
        unique.add(fp)

        if accepted:
            current_structure = cand_structure_relaxed
            current_energy = cand_energy
            energies_accepted.append(current_energy)

            # Save accepted snapshot
            Poscar(current_structure).write_file(
                os.path.join(accepted_dir, f"POSCAR_step_{step:04d}.vasp")
            )

            with open(accepted_log, "a") as f:
                f.write(f"{step}\t{cand_energy:.8f}\t{dE:.8f}\t{accept_prob:.6f}\n")

            print(f"ACCEPTED: E = {cand_energy:.6f} eV, dE = {dE:.6f}, p = {accept_prob:.4f}, u = {u:.4f}")

            # Update best
            if current_energy < best_energy:
                best_energy = current_energy
                best_structure = current_structure
                Poscar(best_structure).write_file(os.path.join(output_dir, "BEST_POSCAR.vasp"))
                print(f"*** New BEST energy: {best_energy:.6f} eV (saved BEST_POSCAR.vasp)")

        else:
            # Save rejected snapshot (unrelaxed or relaxed candidate - we use relaxed one if relax=True)
            Poscar(cand_structure_relaxed).write_file(
                os.path.join(rejected_dir, f"POSCAR_step_{step:04d}.vasp")
            )

            with open(rejected_log, "a") as f:
                f.write(f"{step}\t{cand_energy:.8f}\t{dE:.8f}\t{accept_prob:.6f}\n")

            print(f"REJECTED: E = {cand_energy:.6f} eV, dE = {dE:.6f}, p = {accept_prob:.4f}, u = {u:.4f}")

    # Plot accepted energies trend
    plt.figure()
    plt.plot(range(len(energies_accepted)), energies_accepted, marker="o", linewidth=1)
    plt.xlabel("Accepted move index")
    plt.ylabel("Energy (eV)")
    plt.title("Energy vs Accepted Move (ASE-MLIP MC)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "energy_trend.png"), dpi=200)

    # Also save final best energy
    with open(os.path.join(output_dir, "summary.txt"), "w") as f:
        f.write(f"MC_TEMPERATURE(K): {MC_TEMPERATURE}\n")
        f.write(f"MAX_ITERATIONS: {MAX_ITERATIONS}\n")
        f.write(f"RELAX_BEFORE_ENERGY: {RELAX_BEFORE_ENERGY}\n")
        f.write(f"BestEnergy(eV): {best_energy:.10f}\n")

    print("\nDone.")
    print(f"Best energy: {best_energy:.6f} eV")
    print(f"Outputs saved in: {output_dir}")


# ---------------------------
# Run
# ---------------------------

if __name__ == "__main__":
    monte_carlo_mlip(
        poscar_file="POSCAR",
        scaling_matrix=[3, 3, 1],  # supercell
        substitutions={"Mo": {"Cr": 0.33, "Ti": 0.33}},  # replace Mo sites
        x_elements=["C", "N", "F", "O"],  # protected elements (not substituted; also can be frozen in relax)
        output_dir="mc_output_mlip"
    )
