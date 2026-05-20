#!/usr/bin/env python3
"""
Extract the stress tensor from a VASP vasprun.xml and compute the
hydrostatic pressure and deviatoric stress.

Units
-----
  ase returns stress in eV/Å³ (tension-positive convention).
  Multiply by 160.21766208 to convert to GPa.

Author: Asif Iqbal
Usage:  python3 stress_pressure.py [vasprun.xml] [--example]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import ase.io

_EV_ANG3_TO_GPA = 160.21766208  # 1 eV/Å³ = 160.218 GPa

_RED   = "\033[91m"
_YEL   = "\033[33m"
_RESET = "\033[0m"


# ---------------------------------------------------------------------------
# Pure numerical functions (no I/O — easy to unit-test)
# ---------------------------------------------------------------------------

def hydrostatic_pressure(stress: np.ndarray) -> float:
    """
    Hydrostatic pressure from a 3×3 stress tensor (tension-positive).

    P = −Tr(σ) / 3

    For compressive loading σ_ii < 0 → P > 0.

    Parameters
    ----------
    stress : np.ndarray, shape (3, 3)
        Cauchy stress tensor in eV/Å³ (tension positive).

    Returns
    -------
    float
        P in the same units as *stress*.
    """
    return -float(np.trace(stress)) / 3.0


def deviatoric_stress_matrix(stress: np.ndarray) -> np.ndarray:
    """
    Deviatoric (traceless) part of a 3×3 stress tensor.

    s = σ − (Tr(σ)/3) I  =  σ + P·I

    For isotropic stress σ = σ₀ I this returns the zero matrix.

    Parameters
    ----------
    stress : np.ndarray, shape (3, 3)

    Returns
    -------
    np.ndarray, shape (3, 3)
    """
    p = hydrostatic_pressure(stress)
    return stress + p * np.eye(3)


def deviatoric_stress_voigt(stress_voigt: np.ndarray) -> tuple[np.ndarray, float]:
    """
    Subtract hydrostatic pressure from the diagonal Voigt components.

    Voigt order: [σ₁₁, σ₂₂, σ₃₃, σ₂₃, σ₁₃, σ₁₂]

    P = −mean(σ₁₁, σ₂₂, σ₃₃)
    s_ii = σ_ii + P   (diagonal components only)
    s_ij = σ_ij       (off-diagonal unchanged)

    Parameters
    ----------
    stress_voigt : np.ndarray, shape (6,)

    Returns
    -------
    (dev, p) where dev is the Voigt deviatoric and p is the pressure,
    both in the same units as *stress_voigt*.
    """
    p = -float(np.mean(stress_voigt[:3]))
    dev = stress_voigt.copy()
    dev[:3] += p   # σ_ii + P  (NOT -= )
    return dev, p


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def load_stress(path: Path) -> tuple[np.ndarray, float, np.ndarray]:
    """
    Read a VASP vasprun.xml and return (stress_3x3, volume_Å3, cell_Å).

    The stress is in eV/Å³, tension-positive (ase convention).
    """
    atoms          = ase.io.read(str(path))
    stress_matrix  = atoms.get_stress(voigt=False)
    volume         = float(atoms.get_volume())
    cell           = np.array(atoms.get_cell())
    return stress_matrix, volume, cell


def print_stress_report(path: Path) -> None:
    stress, volume, cell = load_stress(path)

    print(f"\n{'Reading vasprun.xml via ase':*^80}\n")
    print(f"Volume           : {volume:10.4f} Å³")
    atoms = ase.io.read(str(path))
    print(f"{_RED}Lengths & angles : {atoms.get_cell_lengths_and_angles()}{_RESET}")
    print(f"{_RED}Lattice (Å):\n{cell}{_RESET}")
    print("=" * 80)

    print(f"{_YEL}Stress tensor (eV/Å³):\n{stress}{_RESET}")

    p_eVA3 = hydrostatic_pressure(stress)
    p_GPa  = p_eVA3 * _EV_ANG3_TO_GPA
    dev    = deviatoric_stress_matrix(stress)

    print(f"\nHydrostatic pressure  : {p_GPa:10.4f} GPa")
    print(f"Deviatoric stress (eV/Å³):\n{dev}\n")


def manual_kbar_example() -> None:
    """
    Demonstrate manual kbar → GPa conversion with a hard-coded stress block.

    In vasprun.xml stress is in kbar.  VASP stores −σ (compression positive),
    so multiply by −0.1 to get GPa with the tension-positive convention.
    """
    print(f"\n{'Manual kbar→GPa conversion example':_^80}")
    s_kbar = np.array([
        [-0.05224891,  0.00659668,  0.00785320],
        [ 0.00659642, -0.07778282,  0.00069340],
        [ 0.00785295,  0.00069324, -0.07633123],
    ])
    stress_GPa = s_kbar * (-0.1)   # kbar → GPa, flip VASP sign
    print("Stress (GPa):\n", stress_GPa)

    # Voigt ordering: [σ₁₁, σ₂₂, σ₃₃, σ₂₃, σ₁₃, σ₁₂]
    voigt = stress_GPa.reshape(9)[[0, 4, 8, 5, 2, 1]]
    dev, p = deviatoric_stress_voigt(voigt)

    print(f"Hydrostatic pressure  : {p:10.4f} GPa")
    print("Deviatoric stress (Voigt, GPa):\n", dev)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Extract stress tensor and pressure from vasprun.xml."
    )
    p.add_argument(
        "vasprun", nargs="?", default="vasprun.xml",
        help="Path to vasprun.xml (default: vasprun.xml)",
    )
    p.add_argument(
        "--example", action="store_true",
        help="Only run the manual kbar→GPa conversion example",
    )
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    if args.example:
        manual_kbar_example()
    else:
        vasprun = Path(args.vasprun)
        if not vasprun.exists():
            sys.exit(f"ERROR: {vasprun} not found.\n")
        print_stress_report(vasprun)
        manual_kbar_example()
