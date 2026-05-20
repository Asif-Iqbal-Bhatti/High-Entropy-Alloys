# Crystal Deformation — Energy vs Strain

[![License: GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)

Tools for generating Lagrangian-strained crystal structures for elastic constant
calculations in VASP and the exciting DFT code.  The energy-vs-strain approach
fits a polynomial to E(η) to extract elastic constants Cᵢⱼ.

**Author:** Asif Iqbal Bhatti  
**Affiliation:** LSPM CNRS lab Paris France

---

## File overview

```
_strain_core.py              shared numerical kernels (imported by both Python scripts)
lagrangian_strain.py         general deformation — reads VASP CONTCAR, writes POSCARs
cubic_strain.py              cubic-specific deformation — reads exciting input.xml
stress_pressure.py           extract stress tensor & pressure from vasprun.xml
defor_monoclinic_struct.cpp  C++ tool for monoclinic/general deformation along a, b, c
tests/                       pytest suite (76 tests, pure NumPy — no DFT files needed)
```

---

## Dependencies

```
pip install numpy colorama lxml ase
```

For the C++ code: a C++17-capable compiler (`g++ ≥ 7` or `icpc`).

---

## Quick start

### 1. General Lagrangian strain — `lagrangian_strain.py`

Reads a fully relaxed **CONTCAR** (`IBRION=2`, `ISIF=3`) and writes a set of
deformed **POSCAR** files spanning `[−η_max, +η_max]`.

```bash
python3 lagrangian_strain.py [work_dir]
```

The script prompts for:
- Maximum Lagrangian strain η_max (e.g. `0.05` for 5 %)
- Number of strain points (odd number recommended, e.g. `11`)
- Deformation code (see table below)

Output goes into `work_dir/` (default: `workdir/`):

```
workdir/
  POSCAR-01 … POSCAR-11    deformed structures
  strain-01 … strain-11    η value for each step
  INFO-elastic-constants   run parameters
```

**Deformation codes:**

| Code | Voigt pattern | Elastic quantity |
|------|--------------|-----------------|
| 0    | (η, η, η, 0, 0, 0)   | volume strain |
| 1    | (η, 0, 0, 0, 0, 0)   | linear strain along x |
| 2    | (0, η, 0, 0, 0, 0)   | linear strain along y |
| 3    | (0, 0, η, 0, 0, 0)   | linear strain along z |
| 4    | (0, 0, 0, η, 0, 0)   | yz shear |
| 5    | (0, 0, 0, 0, η, 0)   | xz shear |
| 6    | (0, 0, 0, 0, 0, η)   | xy shear |
| 7    | (0, 0, 0, η, η, η)   | shear along ⟨111⟩ |
| 8    | (η, η, 0, 0, 0, 0)   | xy in-plane strain |
| 9    | (η,−η, 0, 0, 0, 0)   | xy in-plane shear |
| 10   | (η, η, η, η, η, η)   | global strain |
| 11–14 | mixed combinations  | mixed strains |

---

### 2. Cubic-specific strain — `cubic_strain.py`

Reads an exciting **input.xml** and writes deformed XML files plus VASP POSCARs
(converted via `ase`).  Implements the three deformation matrices used to
extract B₀, (C₁₁−C₁₂), and C₄₄ for cubic systems.

```bash
python3 cubic_strain.py [work_dir]
```

**Deformation codes (cubic only):**

| Code | Voigt pattern | Elastic quantity |
|------|--------------|-----------------|
| 0 | (η,  η,            η, 0, 0,  0) | 9B₀ = 3C₁₁ + 6C₁₂ |
| 1 | (η, −η, 1/(1−η²)−1, 0, 0,  0) | 2(C₁₁ − C₁₂) |
| 2 | (0,  0, 1/(1−η²)−1, 0, 0, 2η) | 2C₄₄ |

Codes 1 and 2 enforce volume conservation to second order in η.

---

### 3. Stress tensor & pressure — `stress_pressure.py`

Reads a VASP **vasprun.xml** and prints the stress tensor, hydrostatic
pressure, and deviatoric stress.

```bash
python3 stress_pressure.py [vasprun.xml]

# Run only the built-in kbar→GPa conversion example (no file needed):
python3 stress_pressure.py --example
```

**Sign convention (tension-positive, ase standard):**

```
P  = −Tr(σ) / 3
s  = σ + P·I          (deviatoric — traceless)
```

---

### 4. Monoclinic deformation — `defor_monoclinic_struct.cpp`

Generates deformed lattice vectors from a VASP **POSCAR**, varying one axis at
a time or uniformly (volume scaling).  Useful as a quick sanity check or for
non-Lagrangian deformation studies.

**Compile:**
```bash
g++ -std=c++17 -O2 -Wall -o defor defor_monoclinic_struct.cpp
```

**Run:**
```bash
./defor
```

The program prompts for the POSCAR filename, optionally rescales the cell to a
target volume, then asks for the deformation direction (`x`, `y`, `z`, or `v`
for volume), number of steps, and amplitude.  Output is written to
`deformationstruct.dat`.

---

## Running the tests

```bash
pip install pytest
python -m pytest tests/ -v
```

All 76 tests are self-contained (no VASP/exciting files required):

```
tests/test_strain_core.py          core numerical kernels
tests/test_lagrangian_strain.py    POSCAR I/O and deformation codes
tests/test_cubic_strain.py         cubic η matrices, volume conservation
tests/test_stress_pressure.py      pressure and deviatoric stress
```

---

## Theory

The Lagrangian strain tensor η and the linear (infinitesimal) strain ε are
related by:

```
η = ε + ½ ε²
```

This is solved iteratively to recover ε from a specified η.  The deformation
matrix is then:

```
D = I + ε
```

and the deformed lattice vectors are:

```
R' = (D · Rᵀ)ᵀ
```

Strain values should stay in the harmonic regime (η ≲ 5–8 %).  Larger strains
can trigger phase transitions and invalidate the quadratic E(η) fit.

---

## References

[1] exciting code — energy vs strain:
    [http://exciting-code.org/nitrogen-energy-vs-strain-calculations](https://exciting-code.org/)

[2] Materials Project — elasticity calculations:
    https://wiki.materialsproject.org/Elasticity_calculations

[3] ELASTIC code (stress–strain approach):
    https://elastic.readthedocs.io/en/stable/

[4] L.D. Landau, E.M. Lifshitz, *Theory of Elasticity*, Elsevier (1986).
    ISBN: 0750626330

[5] IRelast package — Voigt deformation codes for cubic systems:
    https://doi.org/10.1016/j.jallcom.2017.10.139
