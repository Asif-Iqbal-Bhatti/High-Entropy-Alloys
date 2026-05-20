#!/usr/bin/env python3
"""
Generate Lagrangian-strained POSCAR files from a VASP CONTCAR.

Author: Asif Iqbal
Date:   2019-11-12  (modernised 2026)
Usage:  python3 lagrangian_strain.py [work_dir]
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np
from colorama import Back, Style, init

from _strain_core import (
    build_voigt_strains,
    deform_cell,
    lagrangian_to_linear,
    voigt_to_eta_matrix,
)

init(autoreset=True)

# ---------------------------------------------------------------------------
# Deformation-code catalogue
# ---------------------------------------------------------------------------
DEFORMATION_CODES: dict[int, tuple[str, str]] = {
    # Normal-strain codes: positions 0-2 use 'E'/'e' → Voigt eᵢ = η.
    # Shear-strain codes:  positions 3-5 use '2'     → Voigt eᵢ = 2η,
    #                      so the tensor off-diagonal = eᵢ/2 = η.
    # This matches the ElaStic / exciting convention so that post-processing
    # extracts Cᵢⱼ directly without an extra factor-of-2 correction.
    0:  ("EEE000", "volume strain"),
    1:  ("E00000", "linear strain along x"),
    2:  ("0E0000", "linear strain along y"),
    3:  ("00E000", "linear strain along z"),
    4:  ("000200", "yz shear strain   (tensor η₂₃ = η)"),
    5:  ("000020", "xz shear strain   (tensor η₁₃ = η)"),
    6:  ("000002", "xy shear strain   (tensor η₁₂ = η)"),
    7:  ("000222", "shear along ⟨111⟩ (tensor η₂₃=η₁₃=η₁₂=η)"),
    8:  ("EE0000", "xy in-plane strain"),
    9:  ("Ee0000", "xy in-plane shear strain"),
    10: ("EEE222", "global strain"),
    11: ("E00200", "mixed strain  (e₁, 2e₄)"),
    12: ("E00020", "mixed strain  (e₁, 2e₅)"),
    13: ("E00002", "mixed strain  (e₁, 2e₆)"),
    14: ("EE0200", "mixed strain  (e₁, e₂, 2e₄)"),
}

# ---------------------------------------------------------------------------
# POSCAR I/O
# ---------------------------------------------------------------------------

def read_contcar(path: Path) -> dict:
    """
    Parse a VASP POSCAR/CONTCAR file.

    Returns a dict with keys:
        title, scale, axis_matrix, element_line, atom_number_line,
        nat (list[int]), coord_line, tail_lines (coord type + positions).
    """
    lines = path.read_text().splitlines(keepends=True)
    if len(lines) < 8:
        raise ValueError(f"{path}: fewer than 8 lines; not a valid POSCAR.")

    axis_matrix = np.array(
        [[float(x) for x in lines[i].split()] for i in range(2, 5)],
        dtype=float,
    )
    nat = [int(x) for x in lines[6].split()]

    return {
        "title":           lines[0],
        "scale":           float(lines[1]),
        "axis_matrix":     axis_matrix,
        "element_line":    lines[5],
        "atom_number_line": lines[6],
        "nat":             nat,
        "tail_lines":      lines[7:],   # coord-type line + atomic positions
    }


def write_poscar(
    path: Path,
    title: str,
    scale: float,
    new_axes: np.ndarray,
    element_line: str,
    atom_number_line: str,
    tail_lines: list[str],
) -> None:
    with path.open("w") as f:
        f.write(title)
        f.write(f"{scale:10.8f}\n")
        for row in new_axes:
            f.write(f"{row[0]:22.16f} {row[1]:22.16f} {row[2]:22.16f}\n")
        f.write(element_line)
        f.write(atom_number_line)
        f.writelines(tail_lines)


# ---------------------------------------------------------------------------
# UI helpers
# ---------------------------------------------------------------------------

def print_header() -> None:
    print(f"\n{'Documentation':_^80}")
    print(f"  {'Author':<10}: Asif Iqbal")
    print(f"  {'Date':<10}: 10/02/2021")
    print(f"  {'Usage':<10}: python3 lagrangian_strain.py [work_dir]")
    print("\033[91m  η = ε + ½ε²     (Lagrangian deformation scheme)\033[0m")
    print("\033[91m  r' = (I + ε)·r  (basis-vector transformation)\033[0m")
    print(f"{'END of Documentation':_^80}\n")
    print(Style.RESET_ALL)
    print(Back.GREEN + "-" * 72)
    print(Back.GREEN + "  Deformation codes — Voigt notation")
    print(Back.GREEN + "-" * 72)
    for code, (dc, label) in DEFORMATION_CODES.items():
        color = Back.YELLOW if code in {0, 1, 7} else Back.GREEN
        print(color + f"  {code:>2}  {dc}  {label}")
    print(Back.GREEN + "-" * 72)
    print(Style.RESET_ALL)


def prompt_inputs() -> tuple[float, int, int]:
    maximum_strain = float(input("Enter maximum Lagrangian strain  [e.g. 0.05 for 5%] >>>> "))
    if not (0.0 < maximum_strain <= 1.0):
        sys.exit("ERROR: strain must be in (0, 1].\n")

    strain_points = int(input("Enter number of strain values (odd preferred) >>>> "))
    if not (3 <= strain_points <= 99):
        sys.exit("ERROR: strain_points must be in [3, 99].\n")

    half = strain_points // 2
    print(f"Deformation range: [{-half}, {half}]")

    deformation_code = int(input("\nEnter deformation code >>>> "))
    if deformation_code not in DEFORMATION_CODES:
        sys.exit(f"ERROR: deformation code must be in [0, {max(DEFORMATION_CODES)}].\n")

    return maximum_strain, strain_points, deformation_code


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def elastic_strain(work_directory: str = "workdir") -> None:
    print_header()

    contcar = Path("CONTCAR")
    if not contcar.exists():
        sys.exit("ERROR: CONTCAR not found in the current directory.\n")

    maximum_strain, strain_points, deformation_code = prompt_inputs()
    dc, label = DEFORMATION_CODES[deformation_code]

    poscar = read_contcar(contcar)
    axis_matrix  = poscar["axis_matrix"]
    scale        = poscar["scale"]
    volume0      = abs(np.linalg.det(axis_matrix) * scale**3)

    print(f"Number of atoms : {sum(poscar['nat'])}")
    print(f"Equilibrium vol : {volume0:.6f} Å³")

    work_dir = Path(work_directory)
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir()

    (work_dir / "INFO-elastic-constants").write_text(
        f"Maximum Lagrangian strain       = {maximum_strain}\n"
        f"Number of strain values         = {strain_points}\n"
        f"Volume of equilibrium unit cell = {volume0 * 6.74833:.6f} [a.u.]³\n"
        f"Deformation code                = {deformation_code}\n"
        f"Deformation label               = {dc}  ({label})\n"
    )

    delta    = strain_points - 1
    eta_step = 2.0 * maximum_strain / delta
    half     = strain_points // 2

    hdr = "Vol_D'"
    print(f"\n{'':12s} {'Vol_cell':>12s} {hdr:>14s} {'V/V_D':>14s}")

    for i in range(strain_points):
        eta = i * eta_step - maximum_strain
        (work_dir / f"strain-{i+1:02d}").write_text(f"{eta:10.5f}\n")

        e          = build_voigt_strains(dc, eta)
        eta_mat    = voigt_to_eta_matrix(e)
        eps_mat    = lagrangian_to_linear(eta_mat)
        def_mat    = np.eye(3) + eps_mat
        new_axes   = deform_cell(axis_matrix, def_mat)

        V     = abs(np.linalg.det(new_axes))
        V_def = abs(np.linalg.det(def_mat))
        t, tmp = i + 1, i - half
        print(f"{t:02d}({tmp:3d}) => {V:12.6f} {V_def:14.6f} {V/V_def:14.6f}")

        write_poscar(
            work_dir / f"POSCAR-{t:02d}",
            poscar["title"], scale, new_axes,
            poscar["element_line"], poscar["atom_number_line"], poscar["tail_lines"],
        )

    print()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate Lagrangian-strained POSCAR files from CONTCAR."
    )
    p.add_argument("work_dir", nargs="?", default="workdir",
                   help="Output directory (default: workdir)")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    elastic_strain(args.work_dir)
