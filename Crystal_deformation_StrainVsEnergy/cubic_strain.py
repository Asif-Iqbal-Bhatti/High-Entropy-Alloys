#!/usr/bin/env python3
"""
Generate Lagrangian-strained exciting input.xml / VASP POSCAR files
for cubic crystal systems.

Deformation codes (Voigt notation):
  0 → (η,  η,            η, 0, 0,  0)  volumetric       [9B₀ = 3C₁₁ + 6C₁₂]
  1 → (η, -η,  1/(1-η²)−1, 0, 0,  0)  orthorhombic     [2(C₁₁ − C₁₂)]
  2 → (0,  0,  1/(1-η²)−1, 0, 0, 2η)  monoclinic shear [2C₄₄]

Author: Asif Iqbal
Date:   2020-03-05  (modernised 2026)
Usage:  python3 cubic_strain.py [work_dir]
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from lxml import etree

from _strain_core import (
    build_voigt_strains,
    deform_cell,
    lagrangian_to_linear,
    voigt_to_eta_matrix,
)

_RED   = "\033[91m"
_YEL   = "\033[33m"
_CYAN  = "\033[46m"
_RESET = "\033[0m"

DEFORMATION_CODES: dict[int, tuple[str, str]] = {
    0: ("EEE000", "volumetric strain       [9B₀ = 3C₁₁ + 6C₁₂]"),
    1: ("EeE000", "orthorhombic distortion [2(C₁₁ − C₁₂)]"),
    2: ("00E002", "monoclinic shear        [2C₄₄]"),
}

# ---------------------------------------------------------------------------
# Cubic-specific strain tensor
# ---------------------------------------------------------------------------

def build_eta_matrix_cubic(e: np.ndarray, deformation_code: int) -> np.ndarray:
    """
    Construct the 3×3 Lagrangian strain tensor for a cubic deformation code.

    For codes 1 and 2 the [2,2] element is replaced by 1/(1-e[2]²)-1 to
    enforce volume conservation to second order.

    Parameters
    ----------
    e : np.ndarray, shape (6,)
        Voigt strain components from build_voigt_strains.
    deformation_code : int
        0, 1, or 2.
    """
    eta = np.array([
        [e[0],    e[5]/2, e[4]/2],
        [e[5]/2,  e[1],   e[3]/2],
        [e[4]/2,  e[3]/2, e[2]  ],
    ], dtype=float)

    if deformation_code in (1, 2):
        # Volume-conserving constraint: det(I + η) ≈ 1 to O(η²)
        eta[2, 2] = 1.0 / (1.0 - e[2]**2) - 1.0

    return eta


# ---------------------------------------------------------------------------
# exciting XML I/O
# ---------------------------------------------------------------------------

def read_input_xml(path: Path) -> tuple[etree._ElementTree, np.ndarray, float, list[float]]:
    """
    Parse an exciting input.xml file.

    Returns
    -------
    (doc, axis_matrix, scale, stretch)
        axis_matrix : np.ndarray, shape (3, 3) — rows are basevect
        scale       : float  — crystal/@scale attribute (default 1)
        stretch     : list[float, float, float]  — crystal/@stretch (default [1,1,1])
    """
    doc = etree.parse(str(path))

    scale_attr = doc.xpath("/input/structure/crystal/@scale")
    scale = float(scale_attr[0]) if scale_attr else 1.0

    stretch_attr = doc.xpath("/input/structure/crystal/@stretch")
    stretch = [float(x) for x in stretch_attr[0].split()] if stretch_attr else [1.0, 1.0, 1.0]

    basevects = doc.xpath("//basevect/text()")
    if len(basevects) != 3:
        raise ValueError(f"Expected 3 basevect elements, found {len(basevects)}.")
    axis_matrix = np.array(
        [[float(v) for v in bv.split()] for bv in basevects],
        dtype=float,
    )
    return doc, axis_matrix, scale, stretch


# ---------------------------------------------------------------------------
# UI helpers
# ---------------------------------------------------------------------------

def print_header() -> None:
    print("=" * 80)
    print(f"  {'Author':<10}: Asif Iqbal")
    print(f"  {'Date':<10}: 2020-03-05")
    print(f"  {'Usage':<10}: python3 cubic_strain.py [work_dir]")
    print()
    print("         | η[0]    η[5]/2  η[4]/2 |")
    print("     η = | η[5]/2  η[1]    η[3]/2 |    D' = I + η")
    print("         | η[4]/2  η[3]/2  η[2]   |")
    print("-" * 80)
    for code, (dc, label) in DEFORMATION_CODES.items():
        print(f"{_YEL}  {code}  {dc}  {label}{_RESET}")
    print("=" * 80)


def prompt_inputs() -> tuple[float, int, int]:
    maximum_strain = float(input("Enter maximum Lagrangian strain  [e.g. 0.05 for 5%] >>>> "))
    if not (0.0 < maximum_strain <= 1.0):
        sys.exit("ERROR: strain must be in (0, 1].\n")

    strain_points = int(input("Enter number of strain values (odd preferred) >>>> "))
    if not (3 <= strain_points <= 99):
        sys.exit("ERROR: strain_points must be in [3, 99].\n")

    half = strain_points // 2
    print(f"Deformation range: [{-half}, {half}]")

    deformation_code = int(input("Enter deformation code >>>> "))
    if deformation_code not in DEFORMATION_CODES:
        sys.exit(f"ERROR: deformation code must be one of {list(DEFORMATION_CODES)}.\n")

    return maximum_strain, strain_points, deformation_code


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def cubic_strains(work_directory: str = "workdir") -> None:
    xml_path = Path("input.xml")
    if not xml_path.exists():
        sys.exit("ERROR: input.xml not found in the current directory.\n")

    print_header()
    maximum_strain, strain_points, deformation_code = prompt_inputs()
    dc, label = DEFORMATION_CODES[deformation_code]

    doc, axis_matrix, scale, stretch = read_input_xml(xml_path)
    root = doc.getroot()

    vol0 = abs(
        np.linalg.det(axis_matrix) * scale**3
        * stretch[0] * stretch[1] * stretch[2]
    )

    work_dir = Path(work_directory)
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir()

    (work_dir / "INFO-elastic-constants").write_text(
        f"Maximum Lagrangian strain       = {maximum_strain}\n"
        f"Number of strain values         = {strain_points}\n"
        f"Volume of equilibrium unit cell = {vol0:.6f} [a.u.]³\n"
        f"Deformation code                = {deformation_code}\n"
        f"Deformation label               = {dc}  ({label})\n"
    )

    delta    = strain_points - 1
    eta_step = 2.0 * maximum_strain / delta
    half     = strain_points // 2

    hdr = "Vol_D'"
    print(f"\n{'':12s} {'Vol_cell':>10s} {hdr:>12s} {'V/V0 (%)':>14s}")

    fmt = "%22.16f"

    for i in range(strain_points):
        eta = i * eta_step - maximum_strain
        (work_dir / f"strain-{i+1:02d}").write_text(f"{eta:11.8f}\n")

        e          = build_voigt_strains(dc, eta)
        eta_mat    = build_eta_matrix_cubic(e, deformation_code)
        eps_mat    = lagrangian_to_linear(eta_mat)
        def_mat    = np.eye(3) + eps_mat
        new_axes   = deform_cell(axis_matrix, def_mat)

        V     = abs(np.linalg.det(new_axes))
        V_def = abs(np.linalg.det(def_mat))
        t, tmp = i + 1, i - half
        print(f"{t:02d}({tmp:3d}) => {V:10.6f} {V_def:12.6f} {100*V/vol0:14.2f}")

        # Patch basevect text nodes in-place
        xbv = doc.xpath("//crystal/basevect")
        for j in range(3):
            xbv[j].text = (
                f"{fmt % new_axes[j, 0]}"
                f"{fmt % new_axes[j, 1]}"
                f"{fmt % new_axes[j, 2]} "
            )

        output_xml = work_dir / f"input-{t:02d}.xml"
        output_xml.write_bytes(
            etree.tostring(
                root,
                method="xml",
                pretty_print=True,
                xml_declaration=False,
                encoding="UTF-8",
            )
        )

        subprocess.call([
            "ase", "-T", "convert",
            "-i", "exciting", "-f", str(output_xml),
            "-o", "vasp", str(work_dir / f"POSCAR_{t:02d}"),
            "--write-args", "direct=True", "vasp5=True",
        ])

    print(f"\n{_CYAN}Files converted to VASP POSCAR format by ase.{_RESET}\n")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate cubic Lagrangian-strained exciting XML / POSCAR files."
    )
    p.add_argument("work_dir", nargs="?", default="workdir",
                   help="Output directory (default: workdir)")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    cubic_strains(args.work_dir)
