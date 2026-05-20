#!/usr/bin/env python3
"""
Madelung Potential Calculator — Ewald Summation
================================================
Computes the Madelung (electrostatic) site potential and total energy
for a periodic crystal structure read from a VASP POSCAR/CONTCAR file.

Charges are read from the 4th column of the atomic coordinate block
(in units of elementary charge e, e.g. Na = +1, Cl = -1).

Method
------
Uses pymatgen's EwaldSummation, which decomposes the lattice sum into:
  - Real-space part   : short-ranged, converges exponentially with cutoff
  - Reciprocal part   : long-range, sum over G-vectors
  - Self-energy part  : removes the self-interaction of each Gaussian
  - Point correction  : handles the uniform background for near-neutral cells

This avoids the conditional-convergence problem of direct real-space sums
and gives O(N^{3/2}) scaling with guaranteed accuracy.

Output
------
  - Per-site: charge (e), Madelung potential (V), site energy (eV)
  - Total Madelung energy (eV)
  - Madelung constant (dimensionless, for binary ionic reference)
  - Summary written to madelung_results.tsv

Usage
-----
  python madelung_ewald.py [POSCAR/CONTCAR file]

  If no file is given, looks for CONTCAR then POSCAR in the current directory.

Dependencies
------------
  numpy, pymatgen

Author : Asif Iqbal (refactored with Ewald summation)
"""

import os
import sys
import argparse
import numpy as np

# ---------------------------------------------------------------------------
# Physical constants (CODATA 2018)
# ---------------------------------------------------------------------------
E_CHG  = 1.602176634e-19   # C
EPS0   = 8.8541878128e-12  # F/m
A2M    = 1e-10             # m per Angstrom
KE     = 1.0 / (4.0 * np.pi * EPS0)  # N m^2 / C^2
EV2J   = E_CHG             # 1 eV = 1.602e-19 J


# ===========================================================================
# 1. POSCAR / CONTCAR parser
# ===========================================================================

def read_vasp(filename: str):
    """
    Parse a VASP POSCAR/CONTCAR file.

    Expects charges in the 4th column of the coordinate block (after x y z).
    Returns
    -------
    lattice   : (3, 3) float array, lattice vectors in Angstroms (row vectors)
    frac_coords: (N, 3) float array, fractional coordinates
    cart_coords: (N, 3) float array, Cartesian coordinates in Angstroms
    charges   : (N,)   float array, formal charges in units of e
    species   : list[str], element symbols (or 'X' if not present)
    comment   : str, first line of POSCAR
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Structure file '{filename}' not found.")

    with open(filename, "r") as fh:
        lines = fh.readlines()

    try:
        comment = lines[0].strip()
        scale   = float(lines[1].strip())
        lattice = np.array(
            [line.split() for line in lines[2:5]], dtype=float
        ) * scale

        # Line 5: optional element names (VASP 5+ format)
        # Line 6: atom counts
        # Detect whether line 5 is element names or counts
        line5_tokens = lines[5].split()
        if line5_tokens[0].lstrip("-").replace(".", "", 1).isdigit():
            # Old-style POSCAR: no element names on line 5
            species_names = None
            counts = list(map(int, line5_tokens))
            base   = 6
        else:
            species_names = line5_tokens
            counts = list(map(int, lines[6].split()))
            base   = 7

        total_atoms = sum(counts)

        # Selective dynamics?
        if lines[base].strip().lower().startswith("s"):
            coord_line = base + 1
        else:
            coord_line = base

        coord_type = lines[coord_line].strip().lower()
        coord_line += 1

        coords_raw = []
        charges    = []
        for i in range(total_atoms):
            parts = lines[coord_line + i].split()
            coords_raw.append([float(parts[0]), float(parts[1]), float(parts[2])])
            charges.append(float(parts[3]) if len(parts) >= 4 else 0.0)

        coords_raw = np.array(coords_raw, dtype=float)
        charges    = np.array(charges, dtype=float)

        # Build fractional and Cartesian coordinates
        if coord_type.startswith("d"):          # Direct / fractional
            frac_coords = coords_raw
            cart_coords = frac_coords @ lattice
        else:                                    # Cartesian
            cart_coords = coords_raw * scale
            inv_lat     = np.linalg.inv(lattice)
            frac_coords = cart_coords @ inv_lat

        # Build species list
        if species_names is not None:
            species = []
            for sym, cnt in zip(species_names, counts):
                species.extend([sym] * cnt)
        else:
            species = ["X"] * total_atoms

        return lattice, frac_coords, cart_coords, charges, species, comment

    except Exception as exc:
        raise ValueError(f"Failed to parse VASP file '{filename}': {exc}") from exc


# ===========================================================================
# 2. Ewald summation via pymatgen
# ===========================================================================

def build_pymatgen_structure(lattice, frac_coords, charges, species):
    """
    Construct a pymatgen Structure with oxidation-state-decorated species
    so that EwaldSummation can read formal charges directly.
    """
    try:
        from pymatgen.core import Structure, Lattice, DummySpecies
        from pymatgen.core.periodic_table import Species as PMGSpecies
    except ImportError:
        raise ImportError(
            "pymatgen is required. Install with:  pip install pymatgen"
        )

    pm_lattice = Lattice(lattice)
    pm_species = []
    for sym, q in zip(species, charges):
        try:
            pm_species.append(PMGSpecies(sym, oxidation_state=q))
        except Exception:
            # Fallback to DummySpecies for unknown/custom labels
            pm_species.append(DummySpecies(sym, oxidation_state=q))

    struct = Structure(
        lattice=pm_lattice,
        species=pm_species,
        coords=frac_coords,
        coords_are_cartesian=False,
    )
    return struct


def calculate_ewald(lattice, frac_coords, charges, species, eta=None):
    """
    Compute Madelung site potentials and total energy via Ewald summation.

    Parameters
    ----------
    lattice     : (3, 3) array, lattice vectors in Å
    frac_coords : (N, 3) array, fractional coordinates
    charges     : (N,)   array, formal charges in e
    species     : list[str], element symbols
    eta         : float or None
        Ewald splitting parameter (1/Å²). If None, pymatgen chooses
        an optimal value automatically (~O(N^{1/2}) optimal).

    Returns
    -------
    potentials     : (N,) array, electrostatic potential at each site (Volts)
    site_energies  : (N,) array, electrostatic energy at each site (eV)
    total_energy   : float, total Madelung energy (eV)
    energy_breakdown: dict with real, recip, point, total contributions (eV)
    """
    from pymatgen.analysis.ewald import EwaldSummation

    struct = build_pymatgen_structure(lattice, frac_coords, charges, species)

    kwargs = {}
    if eta is not None:
        kwargs["eta"] = eta

    ew = EwaldSummation(struct, **kwargs)

    # Site energies: E_i = q_i * V_i  =>  V_i = E_i / q_i
    site_energies = np.array([ew.get_site_energy(i) for i in range(len(charges))])

    # Potential in Volts: E_site [eV] / q [e] gives eV/e = Volts
    potentials = np.where(
        np.abs(charges) > 1e-10,
        site_energies / charges,
        0.0,
    )

    total_energy = ew.total_energy  # eV (includes all corrections)

    energy_breakdown = {
        "real_space" : ew.real_space_energy,
        "reciprocal" : ew.reciprocal_space_energy,
        "point_charge": ew.point_energy,
        "total"      : total_energy,
    }

    return potentials, site_energies, total_energy, energy_breakdown


# ===========================================================================
# 3. Madelung constant (for binary ionic reference AB)
# ===========================================================================

def madelung_constant(lattice, frac_coords, charges, species):
    """
    Dimensionless Madelung constant A defined via:
      E_total = -A * N/2 * (q^2 * e^2) / (4 pi eps0 * r0)
    where r0 is the nearest-neighbour distance and N is the number of atoms.

    Only meaningful for a binary compound (two charge types) where the
    formula is well-defined.
    """
    from pymatgen.core import Structure, Lattice

    struct = Structure(Lattice(lattice), species, frac_coords)
    r0 = min(
        struct.get_distance(i, j)
        for i in range(len(charges))
        for j in range(i + 1, len(charges))
        if abs(charges[j]) > 0
    )
    if r0 < 1e-5:
        return None

    n        = len(charges)
    q_ref    = max(abs(charges))
    # Energy unit: q² e² / (4πε₀ r0)  in eV
    # KE in N m² C⁻²; E_CHG in C; A2M in m/Å
    e2_over_4pieps0_r0 = KE * (E_CHG ** 2) / (r0 * A2M) / EV2J  # eV
    energy_unit = q_ref ** 2 * e2_over_4pieps0_r0

    # Use direct real-space sum for a quick estimate of E_total
    # (we already have it from Ewald — caller should pass it in)
    return None  # placeholder; see main() for actual computation


# ===========================================================================
# 4. Output utilities
# ===========================================================================

def print_banner(text: str, width: int = 60):
    print("=" * width)
    print(f"  {text}")
    print("=" * width)


def print_results(charges, species, potentials, site_energies,
                  total_energy, energy_breakdown, r0=None):

    n = len(charges)
    print_banner("Madelung Potential — Ewald Summation Results")

    # ---- per-site table ----
    col = "{:<6} {:<6} {:>12} {:>15} {:>14}"
    print(col.format("Index", "Spec", "Charge (e)", "Potential (V)", "E_site (eV)"))
    print("-" * 60)
    for i in range(n):
        print(col.format(
            i + 1,
            species[i],
            f"{charges[i]:+.4f}",
            f"{potentials[i]:+.6f}",
            f"{site_energies[i]:+.6f}",
        ))

    print("=" * 60)

    # ---- energy breakdown ----
    print(f"\n{'Energy Breakdown':}")
    print(f"  Real-space contribution  : {energy_breakdown['real_space']:+.6f} eV")
    print(f"  Reciprocal contribution  : {energy_breakdown['reciprocal']:+.6f} eV")
    print(f"  Point-charge correction  : {energy_breakdown['point_charge']:+.6f} eV")
    print(f"  Total Madelung energy    : {energy_breakdown['total']:+.6f} eV")

    # ---- Madelung constant (binary case) ----
    unique_charges = np.unique(np.round(charges, 4))
    if len(unique_charges) == 2 and r0 is not None:
        q_ref  = max(abs(unique_charges))
        e2_r0  = KE * (E_CHG ** 2) / (r0 * A2M) / EV2J
        A_mad  = -total_energy / (0.5 * n * q_ref ** 2 * e2_r0)
        print(f"\n  Nearest-neighbour distance : {r0:.4f} Å")
        print(f"  Madelung constant A        : {A_mad:.6f}")
        print(f"  (NaCl reference: 1.747565)")

    print("=" * 60)


def write_tsv(filename, charges, species, potentials, site_energies):
    """Write per-site results to a tab-separated file."""
    with open(filename, "w") as fh:
        fh.write("index\tspecies\tcharge_e\tpotential_V\tsite_energy_eV\n")
        for i, (sym, q, v, e) in enumerate(
            zip(species, charges, potentials, site_energies), start=1
        ):
            fh.write(f"{i}\t{sym}\t{q:+.6f}\t{v:+.6f}\t{e:+.6f}\n")
    print(f"\nResults written to '{filename}'")


# ===========================================================================
# 5. Validation: compare with textbook Madelung constant
# ===========================================================================

def nacl_validation_test():
    """
    Build a 2-atom NaCl primitive cell and check that the Madelung constant
    is ~1.7476. Returns True if the test passes within 1%.
    """
    a  = 5.6402  # Å  (conventional cubic cell parameter)
    a2 = a / 2.0
    lattice = np.array([
        [0.0, a2,  a2 ],
        [a2,  0.0, a2 ],
        [a2,  a2,  0.0],
    ])
    frac_coords = np.array([[0.0, 0.0, 0.0],   # Na
                             [0.5, 0.5, 0.5]])  # Cl
    charges     = np.array([+1.0, -1.0])
    species     = ["Na", "Cl"]

    potentials, site_energies, total_energy, _ = calculate_ewald(
        lattice, frac_coords, charges, species
    )

    r0     = a2                               # NN distance in FCC NaCl
    q_ref  = 1.0
    e2_r0  = KE * (E_CHG ** 2) / (r0 * A2M) / EV2J
    n      = 2
    A_mad  = -total_energy / (0.5 * n * q_ref ** 2 * e2_r0)
    ref    = 1.7475645946                     # Madelung constant for NaCl

    ok = abs(A_mad - ref) / ref < 0.01
    print(f"[Validation] NaCl Madelung constant: {A_mad:.6f}  "
          f"(ref: {ref:.6f})  {'PASS ✓' if ok else 'FAIL ✗'}")
    return ok


# ===========================================================================
# 6. Main entry point
# ===========================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Madelung potential via Ewald summation from VASP POSCAR/CONTCAR."
    )
    p.add_argument(
        "poscar", nargs="?", default=None,
        help="Path to POSCAR or CONTCAR (default: auto-detect in CWD)",
    )
    p.add_argument(
        "--eta", type=float, default=None,
        help="Ewald splitting parameter η (1/Å²). Default: auto-optimised.",
    )
    p.add_argument(
        "--output", default="madelung_results.tsv",
        help="TSV output filename (default: madelung_results.tsv)",
    )
    p.add_argument(
        "--validate", action="store_true",
        help="Run NaCl validation test before processing.",
    )
    p.add_argument(
        "--no-tsv", action="store_true",
        help="Skip writing TSV output.",
    )
    return p.parse_args()


def find_input_file(given: str | None) -> str:
    if given is not None:
        return given
    for candidate in ("CONTCAR", "POSCAR"):
        if os.path.exists(candidate):
            return candidate
    raise FileNotFoundError(
        "No POSCAR/CONTCAR found in the current directory. "
        "Pass a filename explicitly."
    )


def main():
    args = parse_args()

    # Optional NaCl validation
    if args.validate:
        print_banner("NaCl Validation Test")
        passed = nacl_validation_test()
        if not passed:
            print("WARNING: Validation failed — check your pymatgen installation.")
        print()

    # Locate and parse input
    input_file = find_input_file(args.poscar)
    print(f"Reading: {input_file}")

    lattice, frac_coords, cart_coords, charges, species, comment = read_vasp(input_file)
    n = len(charges)
    print(f"Comment : {comment}")
    print(f"Atoms   : {n}")
    print(f"Net charge: {np.sum(charges):.4f} e")

    net = np.sum(charges)
    if abs(net) > 1e-3:
        print(
            f"WARNING: System net charge = {net:.4f} e. "
            "EwaldSummation adds a uniform compensating background automatically."
        )

    # Nearest-neighbour distance (for Madelung constant)
    try:
        from pymatgen.core import Structure, Lattice
        struct_tmp = Structure(Lattice(lattice), species, frac_coords)
        r0 = min(
            struct_tmp.get_distance(i, j)
            for i in range(n)
            for j in range(i + 1, n)
        )
    except Exception:
        r0 = None

    # Run Ewald
    print("\nRunning Ewald summation...")
    potentials, site_energies, total_energy, breakdown = calculate_ewald(
        lattice, frac_coords, charges, species, eta=args.eta
    )

    # Print results
    print()
    print_results(charges, species, potentials, site_energies,
                  total_energy, breakdown, r0=r0)

    # Save TSV
    if not args.no_tsv:
        write_tsv(args.output, charges, species, potentials, site_energies)


if __name__ == "__main__":
    main()
