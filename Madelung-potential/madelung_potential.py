#!/usr/bin/env python3
"""
Madelung Potential Calculator — Ewald Summation
================================================
Computes the electrostatic (Madelung) site potential and total energy
for a periodic crystal structure read from a VASP POSCAR/CONTCAR file.

Charge sources (choose one via --charges flag)
-----------------------------------------------
  poscar   : 4th column of the POSCAR/CONTCAR coordinate block.
             Charges must be in units of elementary charge e
             (e.g. Na = +1, Cl = -1).  [default]

  bader    : Bader charges from a VASP Bader analysis run.
             Requires ACF.dat (or --bader-file) and POTCAR (or --zvals)
             in the working directory.
             Ionic charge q_i = Z_i - Q_Bader_i
             where Z_i is the ZVAL from POTCAR and Q_Bader_i is the
             integrated Bader electron count from ACF.dat.

  mulliken : Mulliken or NBO partial charges from an external file.
             One charge per line, same atom ordering as the POSCAR.
             Use --charges-file to specify the file path.

Method
------
Uses pymatgen's EwaldSummation, which decomposes the lattice sum into:
  - Real-space part   : short-ranged, converges exponentially with cutoff
  - Reciprocal part   : long-range, sum over G-vectors
  - Self-energy part  : removes the self-interaction of each Gaussian
  - Point correction  : handles the uniform background for near-neutral cells

Output
------
  - Per-site: species, charge (e), Madelung potential (V), site energy (eV)
  - Total Madelung energy (eV) with real/reciprocal/point breakdown
  - Madelung constant A (for binary compounds)
  - TSV file: madelung_results.tsv  (disable with --no-tsv)

Usage examples
--------------
  # Formal charges from POSCAR column 4 (default)
  python madelung_ewald.py POSCAR

  # Bader charges — auto-reads ACF.dat + POTCAR in CWD
  python madelung_ewald.py POSCAR --charges bader

  # Bader charges — explicit files
  python madelung_ewald.py POSCAR --charges bader \
      --bader-file /path/to/ACF.dat --zvals 11 17

  # Mulliken/NBO charges from a plain text file (one value per line)
  python madelung_ewald.py POSCAR --charges mulliken \
      --charges-file nbo_charges.txt

  # Run NaCl validation and use custom Ewald eta
  python madelung_ewald.py POSCAR --validate --eta 0.5

Dependencies
------------
  numpy, pymatgen

Author : Asif Iqbal
"""

import os
import re
import sys
import argparse
import numpy as np

# ---------------------------------------------------------------------------
# Physical constants (CODATA 2018)
# ---------------------------------------------------------------------------
E_CHG = 1.602176634e-19    # C  (elementary charge)
EPS0  = 8.8541878128e-12   # F/m
A2M   = 1e-10              # m per Ångström
KE    = 1.0 / (4.0 * np.pi * EPS0)   # N m² C⁻²
EV2J  = E_CHG              # 1 eV in Joules

# Known nuclear valence charges (ZVAL) for common VASP PAW potentials.
# Used as a fallback when POTCAR is absent and --zvals is not supplied.
# Values match the standard PAW_PBE potentials shipped with VASP 6.
ZVAL_TABLE = {
    "H":  1,  "He": 2,  "Li": 1,  "Be": 2,  "B":  3,  "C":  4,
    "N":  5,  "O":  6,  "F":  7,  "Ne": 8,  "Na": 1,  "Mg": 2,
    "Al": 3,  "Si": 4,  "P":  5,  "S":  6,  "Cl": 7,  "Ar": 8,
    "K":  1,  "Ca": 2,  "Sc": 3,  "Ti": 4,  "V":  5,  "Cr": 6,
    "Mn": 7,  "Fe": 8,  "Co": 9,  "Ni": 10, "Cu": 11, "Zn": 12,
    "Ga": 3,  "Ge": 4,  "As": 5,  "Se": 6,  "Br": 7,  "Kr": 8,
    "Rb": 1,  "Sr": 2,  "Y":  3,  "Zr": 4,  "Nb": 5,  "Mo": 6,
    "Tc": 7,  "Ru": 8,  "Rh": 9,  "Pd": 10, "Ag": 11, "Cd": 12,
    "In": 3,  "Sn": 4,  "Sb": 5,  "Te": 6,  "I":  7,  "Xe": 8,
    "Cs": 1,  "Ba": 2,  "La": 11, "Ce": 12, "Pr": 13, "Nd": 14,
    "Pm": 15, "Sm": 16, "Eu": 17, "Gd": 18, "Tb": 19, "Dy": 20,
    "Ho": 21, "Er": 22, "Tm": 23, "Yb": 24, "Lu": 25, "Hf": 4,
    "Ta": 5,  "W":  6,  "Re": 7,  "Os": 8,  "Ir": 9,  "Pt": 10,
    "Au": 11, "Hg": 12, "Tl": 3,  "Pb": 4,  "Bi": 5,  "Po": 6,
}


# ===========================================================================
# Section 1 — POSCAR / CONTCAR parser
# ===========================================================================

def read_vasp(filename: str):
    """
    Parse a VASP POSCAR/CONTCAR file.

    Charges from column 4 are read when present; otherwise set to 0.
    Handles both VASP 4 (no element line) and VASP 5+ (element line present),
    and optional Selective Dynamics blocks.

    Returns
    -------
    lattice     : (3,3) float array  — lattice vectors in Å (row vectors)
    frac_coords : (N,3) float array  — fractional coordinates
    cart_coords : (N,3) float array  — Cartesian coordinates in Å
    charges     : (N,)  float array  — column-4 charges in e (0 if absent)
    species     : list[str]          — element symbols per atom
    counts      : list[int]          — atom counts per species
    comment     : str                — POSCAR title line
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: '{filename}'")

    with open(filename) as fh:
        lines = fh.readlines()

    try:
        comment = lines[0].strip()
        scale   = float(lines[1].strip())
        lattice = np.array([l.split() for l in lines[2:5]], dtype=float) * scale

        # Distinguish VASP 4 (no element names) from VASP 5+
        tok5 = lines[5].split()
        if tok5[0].lstrip("-").replace(".", "", 1).isdigit():
            species_names = None
            counts = list(map(int, tok5))
            base   = 6
        else:
            species_names = tok5
            counts = list(map(int, lines[6].split()))
            base   = 7

        n_atoms = sum(counts)

        # Optional Selective Dynamics line
        coord_line = base + 1 if lines[base].strip().lower().startswith("s") else base
        coord_type = lines[coord_line].strip().lower()
        coord_line += 1

        coords_raw, charges = [], []
        for i in range(n_atoms):
            parts = lines[coord_line + i].split()
            coords_raw.append(parts[:3])
            charges.append(float(parts[3]) if len(parts) >= 4 else 0.0)

        coords_raw = np.array(coords_raw, dtype=float)
        charges    = np.array(charges,    dtype=float)

        if coord_type.startswith("d"):
            frac_coords = coords_raw
            cart_coords = frac_coords @ lattice
        else:
            cart_coords = coords_raw * scale
            frac_coords = cart_coords @ np.linalg.inv(lattice)

        # Build per-atom species list
        if species_names is not None:
            species = []
            for sym, cnt in zip(species_names, counts):
                species.extend([sym] * cnt)
        else:
            species = ["X"] * n_atoms

        return lattice, frac_coords, cart_coords, charges, species, counts, comment

    except Exception as exc:
        raise ValueError(f"Failed to parse '{filename}': {exc}") from exc


# ===========================================================================
# Section 2 — Bader charge reader
# ===========================================================================

def read_potcar_zvals(potcar_path: str) -> list[float]:
    """
    Extract ZVAL (valence electron count) for each species block in a POTCAR.
    Returns a list in the same order as the species in the POTCAR.
    """
    zvals = []
    with open(potcar_path) as fh:
        for line in fh:
            if "ZVAL" in line:
                # e.g.   ZVAL =   1.000
                m = re.search(r"ZVAL\s*=\s*([\d.]+)", line)
                if m:
                    zvals.append(float(m.group(1)))
    if not zvals:
        raise ValueError(
            f"No ZVAL entries found in '{potcar_path}'. "
            "Check that it is a valid POTCAR file."
        )
    return zvals


def build_zvals_array(species: list[str],
                      counts: list[int],
                      potcar_path: str | None,
                      explicit_zvals: list[float] | None) -> np.ndarray:
    """
    Assemble a per-atom ZVAL array from one of three sources (priority order):
      1. explicit_zvals  — provided via --zvals on the command line
      2. POTCAR          — parsed from the POTCAR file
      3. ZVAL_TABLE      — built-in fallback table (common PAW_PBE potentials)

    Parameters
    ----------
    species        : per-atom element symbols
    counts         : atom counts per unique species (same order as POTCAR blocks)
    potcar_path    : path to POTCAR, or None
    explicit_zvals : list of ZVAL values (one per unique species), or None

    Returns
    -------
    zvals_per_atom : (N,) float array
    """
    unique_species = []
    for sym, cnt in zip(
        # reconstruct unique species list from counts + species
        _unique_ordered(species), counts
    ):
        unique_species.append(sym)

    if explicit_zvals is not None:
        if len(explicit_zvals) != len(unique_species):
            raise ValueError(
                f"--zvals has {len(explicit_zvals)} values but structure has "
                f"{len(unique_species)} unique species: {unique_species}"
            )
        zval_map = dict(zip(unique_species, explicit_zvals))
        source = "command line"

    elif potcar_path is not None and os.path.exists(potcar_path):
        raw = read_potcar_zvals(potcar_path)
        if len(raw) != len(unique_species):
            raise ValueError(
                f"POTCAR has {len(raw)} species blocks but structure has "
                f"{len(unique_species)}: {unique_species}. "
                "Check that POTCAR and POSCAR species order matches."
            )
        zval_map = dict(zip(unique_species, raw))
        source = f"POTCAR ({potcar_path})"

    else:
        missing = [s for s in unique_species if s not in ZVAL_TABLE]
        if missing:
            raise ValueError(
                f"Species {missing} not in the built-in ZVAL table and no "
                "POTCAR found. Supply --zvals or --potcar explicitly."
            )
        zval_map = {s: ZVAL_TABLE[s] for s in unique_species}
        source = "built-in ZVAL table"

    print(f"  ZVAL source  : {source}")
    for sym in unique_species:
        print(f"    {sym:4s}  ZVAL = {zval_map[sym]:.1f}")

    return np.array([zval_map[s] for s in species], dtype=float)


def _unique_ordered(seq):
    """Return unique elements of seq preserving first-occurrence order."""
    seen, out = set(), []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def read_acf(acf_path: str, n_atoms: int) -> np.ndarray:
    """
    Parse the ACF.dat file produced by the Henkelman Bader analysis code.

    Format (columns): index  X  Y  Z  CHARGE  MIN_DIST  ATOMIC_VOL
    Returns the CHARGE column (Bader electron counts) as a float array.
    """
    if not os.path.exists(acf_path):
        raise FileNotFoundError(
            f"Bader output file '{acf_path}' not found. "
            "Run the Bader analysis code first (https://theory.cm.utexas.edu/henkelman/code/bader/)."
        )

    bader_charges = []
    with open(acf_path) as fh:
        for line in fh:
            stripped = line.strip()
            # Skip header/footer/separator lines
            if (not stripped
                    or stripped.startswith("#")
                    or stripped.startswith("-")
                    or stripped.upper().startswith("VACUUM")
                    or stripped.upper().startswith("NUMBER")):
                continue
            parts = stripped.split()
            # Data lines: first token is integer atom index
            if not parts[0].isdigit():
                continue
            try:
                bader_charges.append(float(parts[4]))   # column 5 = CHARGE
            except (IndexError, ValueError):
                continue

    if len(bader_charges) != n_atoms:
        raise ValueError(
            f"ACF.dat contains {len(bader_charges)} atom entries but "
            f"POSCAR has {n_atoms} atoms. Check that they match."
        )

    return np.array(bader_charges, dtype=float)


def get_bader_ionic_charges(species: list[str],
                             counts: list[int],
                             acf_path: str,
                             potcar_path: str | None,
                             explicit_zvals: list[float] | None) -> np.ndarray:
    """
    Compute ionic charges from Bader analysis:
        q_i = ZVAL_i - Q_Bader_i

    Parameters
    ----------
    species        : per-atom element symbols (length N)
    counts         : atom counts per unique species
    acf_path       : path to ACF.dat
    potcar_path    : path to POTCAR (None to fall back to ZVAL_TABLE)
    explicit_zvals : per-species ZVAL list from CLI, or None

    Returns
    -------
    ionic_charges : (N,) float array in units of e
    """
    n = len(species)
    print(f"\nReading Bader charges from '{acf_path}'")
    bader_electrons = read_acf(acf_path, n)
    zvals           = build_zvals_array(species, counts, potcar_path, explicit_zvals)
    ionic_charges   = zvals - bader_electrons

    print(f"  Total Bader electrons : {bader_electrons.sum():.4f}")
    print(f"  Total ZVAL electrons  : {zvals.sum():.4f}")
    print(f"  Net ionic charge      : {ionic_charges.sum():+.4f} e")

    return ionic_charges


# ===========================================================================
# Section 3 — External partial charge reader (Mulliken / NBO / custom)
# ===========================================================================

def read_external_charges(filepath: str, n_atoms: int) -> np.ndarray:
    """
    Read partial charges from a plain-text file with one numeric value per line.
    Comment lines starting with '#' are ignored.
    The file must contain exactly n_atoms charge values.

    Compatible with Mulliken, NBO, Hirshfeld, CM5, or any other scheme
    where charges have already been computed externally.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Charges file not found: '{filepath}'")

    charges = []
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                charges.append(float(line.split()[0]))
            except ValueError:
                continue

    if len(charges) != n_atoms:
        raise ValueError(
            f"'{filepath}' contains {len(charges)} values but structure "
            f"has {n_atoms} atoms."
        )

    return np.array(charges, dtype=float)


# ===========================================================================
# Section 4 — Ewald summation via pymatgen
# ===========================================================================

def build_pymatgen_structure(lattice, frac_coords, charges, species):
    """
    Build a pymatgen Structure with oxidation-state-decorated species.
    DummySpecies is used for unrecognised element symbols.
    """
    try:
        from pymatgen.core import Structure, Lattice, DummySpecies
        from pymatgen.core.periodic_table import Species as PMGSpecies
    except ImportError:
        raise ImportError("Install pymatgen:  pip install pymatgen")

    pm_lattice = Lattice(lattice)
    pm_species = []
    for sym, q in zip(species, charges):
        try:
            pm_species.append(PMGSpecies(sym, oxidation_state=q))
        except Exception:
            pm_species.append(DummySpecies(sym, oxidation_state=q))

    return Structure(
        lattice=pm_lattice,
        species=pm_species,
        coords=frac_coords,
        coords_are_cartesian=False,
    )


def calculate_ewald(lattice, frac_coords, charges, species, eta=None):
    """
    Ewald summation for electrostatic site potentials and total energy.

    Parameters
    ----------
    lattice     : (3,3) array  — lattice vectors in Å
    frac_coords : (N,3) array  — fractional coordinates
    charges     : (N,)  array  — charges in e
    species     : list[str]    — element symbols
    eta         : float|None   — Ewald splitting parameter (1/Å²);
                                 None → pymatgen auto-optimises

    Returns
    -------
    potentials      : (N,) array  — electrostatic potential at each site (V)
    site_energies   : (N,) array  — electrostatic site energy (eV)
    total_energy    : float       — total Madelung energy (eV)
    energy_breakdown: dict        — real/reciprocal/point/total (eV)
    """
    from pymatgen.analysis.ewald import EwaldSummation

    struct = build_pymatgen_structure(lattice, frac_coords, charges, species)
    kwargs = {} if eta is None else {"eta": eta}
    ew     = EwaldSummation(struct, **kwargs)

    site_energies = np.array([ew.get_site_energy(i) for i in range(len(charges))])

    # V_i [Volt] = E_i [eV] / q_i [e]   (eV/e = V)
    potentials = np.where(
        np.abs(charges) > 1e-10,
        site_energies / charges,
        0.0,
    )

    total_energy     = ew.total_energy
    energy_breakdown = {
        "real_space" : ew.real_space_energy,
        "reciprocal" : ew.reciprocal_space_energy,
        "point_charge": ew.point_energy,
        "total"      : total_energy,
    }
    return potentials, site_energies, total_energy, energy_breakdown


# ===========================================================================
# Section 5 — Output utilities
# ===========================================================================

def nearest_neighbour(lattice, frac_coords, species):
    """Return the shortest inter-atomic distance in Å."""
    try:
        from pymatgen.core import Structure, Lattice
        st = Structure(Lattice(lattice), species, frac_coords)
        n  = len(species)
        return min(
            st.get_distance(i, j)
            for i in range(n) for j in range(i + 1, n)
        )
    except Exception:
        return None


def madelung_constant_binary(total_energy_ev, charges, r0):
    """
    Compute the dimensionless Madelung constant for a binary compound:
        E_total = -A * (N/2) * (q² e²) / (4πε₀ r₀)
    Only valid when there are exactly two unique charge magnitudes.
    """
    unique_q = np.unique(np.abs(np.round(charges, 4)))
    if len(unique_q) != 2 or r0 is None or r0 < 1e-5:
        return None
    n     = len(charges)
    q_ref = max(unique_q)
    e2_r0 = KE * (E_CHG ** 2) / (r0 * A2M) / EV2J   # eV
    A     = -total_energy_ev / (0.5 * n * q_ref ** 2 * e2_r0)
    return A


def print_banner(text, width=62):
    print("=" * width)
    print(f"  {text}")
    print("=" * width)


def print_results(charges, species, potentials, site_energies,
                  total_energy, energy_breakdown, charge_label,
                  r0=None):

    print_banner(f"Electrostatic Potential — Ewald Summation  [{charge_label}]")

    hdr = "{:<6} {:<6} {:>13} {:>16} {:>15}"
    row = "{:<6} {:<6} {:>+13.5f} {:>+16.6f} {:>+15.6f}"
    print(hdr.format("Index", "Spec", "Charge (e)", "Potential (V)", "E_site (eV)"))
    print("-" * 62)
    for i, (sym, q, v, e) in enumerate(
        zip(species, charges, potentials, site_energies), start=1
    ):
        print(row.format(i, sym, q, v, e))
    print("=" * 62)

    print(f"\nEnergy breakdown")
    print(f"  Real-space part      : {energy_breakdown['real_space']:+.6f} eV")
    print(f"  Reciprocal-space part: {energy_breakdown['reciprocal']:+.6f} eV")
    print(f"  Point-charge term    : {energy_breakdown['point_charge']:+.6f} eV")
    print(f"  Total Madelung energy: {energy_breakdown['total']:+.6f} eV")

    A = madelung_constant_binary(total_energy, charges, r0)
    if A is not None:
        print(f"\n  Nearest-neighbour r₀ : {r0:.4f} Å")
        print(f"  Madelung constant A  : {A:.6f}")
        print(f"  (NaCl reference      : 1.747565)")

    print("=" * 62)


def write_tsv(path, charges, species, potentials, site_energies, charge_label):
    header = (
        f"# Charge source: {charge_label}\n"
        "index\tspecies\tcharge_e\tpotential_V\tsite_energy_eV\n"
    )
    with open(path, "w") as fh:
        fh.write(header)
        for i, (sym, q, v, e) in enumerate(
            zip(species, charges, potentials, site_energies), start=1
        ):
            fh.write(f"{i}\t{sym}\t{q:+.6f}\t{v:+.6f}\t{e:+.6f}\n")
    print(f"\nResults written to '{path}'")


# ===========================================================================
# Section 6 — NaCl validation test
# ===========================================================================

def nacl_validation_test():
    """
    Build a 2-atom NaCl primitive cell and verify that the Madelung constant
    reproduces the textbook value 1.7475645946 to within 0.01 %.
    """
    a2 = 5.6402 / 2.0
    lattice     = np.array([[0, a2, a2], [a2, 0, a2], [a2, a2, 0]], dtype=float)
    frac_coords = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
    charges     = np.array([+1.0, -1.0])
    species     = ["Na", "Cl"]

    _, _, total_energy, _ = calculate_ewald(lattice, frac_coords, charges, species)

    r0    = a2
    e2_r0 = KE * (E_CHG ** 2) / (r0 * A2M) / EV2J
    A     = -total_energy / (0.5 * 2 * 1.0 ** 2 * e2_r0)
    ref   = 1.7475645946

    ok = abs(A - ref) / ref < 1e-4
    print(f"[Validation] NaCl Madelung constant : {A:.6f}  "
          f"(ref: {ref:.6f})  {'PASS ✓' if ok else 'FAIL ✗'}")
    return ok


# ===========================================================================
# Section 7 — CLI and main
# ===========================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Madelung (electrostatic) potential via Ewald summation.\n"
            "Supports formal, Bader, Mulliken, and NBO charges."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    p.add_argument(
        "poscar", nargs="?", default=None,
        help="POSCAR or CONTCAR path. Auto-detects CONTCAR/POSCAR in CWD if omitted.",
    )

    # ---- charge source ----
    cg = p.add_argument_group("charge source (choose one)")
    cg.add_argument(
        "--charges",
        choices=["poscar", "bader", "mulliken"],
        default="poscar",
        metavar="MODE",
        help=(
            "Charge source: 'poscar' (column 4, default), "
            "'bader' (ACF.dat + POTCAR/ZVAL_TABLE), "
            "'mulliken' (external file, also works for NBO/Hirshfeld)."
        ),
    )
    cg.add_argument(
        "--bader-file", default="ACF.dat",
        metavar="PATH",
        help="Path to ACF.dat from Bader analysis (default: ACF.dat in CWD).",
    )
    cg.add_argument(
        "--potcar", default="POTCAR",
        metavar="PATH",
        help="Path to POTCAR for ZVAL extraction (default: POTCAR in CWD).",
    )
    cg.add_argument(
        "--zvals", nargs="+", type=float, default=None,
        metavar="Z",
        help=(
            "Manual ZVAL per unique species in POSCAR order, e.g. "
            "--zvals 11 17 for Na Cl. Overrides POTCAR and built-in table."
        ),
    )
    cg.add_argument(
        "--charges-file", default=None,
        metavar="PATH",
        help="Path to external charges file for --charges mulliken (one value per line).",
    )

    # ---- Ewald options ----
    eg = p.add_argument_group("Ewald options")
    eg.add_argument(
        "--eta", type=float, default=None,
        help="Ewald splitting parameter η in 1/Å². Default: auto-optimised.",
    )

    # ---- output ----
    og = p.add_argument_group("output")
    og.add_argument(
        "--output", default="madelung_results.tsv",
        help="TSV output filename (default: madelung_results.tsv).",
    )
    og.add_argument(
        "--no-tsv", action="store_true",
        help="Skip TSV output.",
    )
    og.add_argument(
        "--validate", action="store_true",
        help="Run NaCl Madelung constant validation before processing.",
    )

    return p.parse_args()


def find_poscar(given):
    if given:
        return given
    for f in ("CONTCAR", "POSCAR"):
        if os.path.exists(f):
            return f
    raise FileNotFoundError(
        "No POSCAR/CONTCAR in CWD. Pass a filename explicitly."
    )


def main():
    args = parse_args()

    # --- optional validation ---
    if args.validate:
        print_banner("NaCl Validation Test")
        ok = nacl_validation_test()
        if not ok:
            print("WARNING: Validation failed — check your pymatgen installation.")
        print()

    # --- parse structure ---
    poscar = find_poscar(args.poscar)
    print(f"Reading structure : {poscar}")
    lattice, frac_coords, cart_coords, poscar_charges, species, counts, comment = \
        read_vasp(poscar)
    n = len(species)
    print(f"Comment           : {comment}")
    print(f"Atoms             : {n}  ({', '.join(f'{s}×{c}' for s,c in zip(_unique_ordered(species), counts))})")

    # --- resolve charges ---
    mode = args.charges.lower()

    if mode == "poscar":
        charges      = poscar_charges
        charge_label = "formal/POSCAR column-4"
        if np.all(charges == 0.0):
            print(
                "\nWARNING: All POSCAR charges are zero. "
                "Add formal charges as the 4th coordinate column,\n"
                "or use --charges bader / --charges mulliken."
            )

    elif mode == "bader":
        potcar_path = args.potcar if os.path.exists(args.potcar) else None
        charges = get_bader_ionic_charges(
            species, counts,
            acf_path       = args.bader_file,
            potcar_path    = potcar_path,
            explicit_zvals = args.zvals,
        )
        charge_label = f"Bader ionic (ACF: {args.bader_file})"

    elif mode == "mulliken":
        if args.charges_file is None:
            raise ValueError(
                "--charges mulliken requires --charges-file <path>."
            )
        charges = read_external_charges(args.charges_file, n)
        charge_label = f"Mulliken/NBO ({args.charges_file})"
        print(f"\nExternal charges read from '{args.charges_file}'")

    else:
        raise ValueError(f"Unknown charge mode: {mode}")

    # --- neutrality check ---
    net = charges.sum()
    print(f"\nNet charge        : {net:+.4f} e")
    if abs(net) > 0.05:
        print(
            "NOTE: System is not neutral. EwaldSummation applies a uniform\n"
            "      compensating background (jellium) automatically."
        )

    # --- nearest-neighbour distance ---
    r0 = nearest_neighbour(lattice, frac_coords, species)

    # --- run Ewald ---
    print("\nRunning Ewald summation ...")
    potentials, site_energies, total_energy, breakdown = calculate_ewald(
        lattice, frac_coords, charges, species, eta=args.eta
    )

    # --- print and save ---
    print()
    print_results(charges, species, potentials, site_energies,
                  total_energy, breakdown, charge_label, r0=r0)

    if not args.no_tsv:
        write_tsv(args.output, charges, species,
                  potentials, site_energies, charge_label)


if __name__ == "__main__":
    main()
