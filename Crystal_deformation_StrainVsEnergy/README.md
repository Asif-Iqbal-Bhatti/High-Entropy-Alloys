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
tests/                       pytest suite (80 tests, pure NumPy — no DFT files needed)
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

| Code | Voigt pattern | Tensor shear | Elastic quantity |
|------|--------------|-------------|-----------------|
| 0    | (η,  η,  η,  0,  0,  0)  | —      | volume strain |
| 1    | (η,  0,  0,  0,  0,  0)  | —      | linear strain along x |
| 2    | (0,  η,  0,  0,  0,  0)  | —      | linear strain along y |
| 3    | (0,  0,  η,  0,  0,  0)  | —      | linear strain along z |
| 4    | (0,  0,  0, 2η,  0,  0)  | η₂₃=η | yz shear |
| 5    | (0,  0,  0,  0, 2η,  0)  | η₁₃=η | xz shear |
| 6    | (0,  0,  0,  0,  0, 2η)  | η₁₂=η | xy shear |
| 7    | (0,  0,  0, 2η, 2η, 2η)  | all=η  | shear along ⟨111⟩ |
| 8    | (η,  η,  0,  0,  0,  0)  | —      | xy in-plane strain |
| 9    | (η, −η,  0,  0,  0,  0)  | —      | xy in-plane shear |
| 10   | (η,  η,  η, 2η, 2η, 2η)  | all=η  | global strain |
| 11–14 | mixed normal + shear    | —      | mixed strains |

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

---

### 3. Stress tensor & pressure — `stress_pressure.py`

Reads a VASP **vasprun.xml** and prints the stress tensor, hydrostatic
pressure, and deviatoric stress.

```bash
python3 stress_pressure.py [vasprun.xml]

# Run only the built-in kbar→GPa conversion example (no file needed):
python3 stress_pressure.py --example
```

---

### 4. Monoclinic deformation — `defor_monoclinic_struct.cpp`

Generates deformed lattice vectors from a VASP **POSCAR**, varying one axis at
a time or uniformly (volume scaling).

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

All 80 tests are self-contained (no VASP/exciting files required):

```
tests/test_strain_core.py          core numerical kernels
tests/test_lagrangian_strain.py    POSCAR I/O and deformation codes
tests/test_cubic_strain.py         cubic η matrices and volume conservation
tests/test_stress_pressure.py      pressure and deviatoric stress
```

---

## Theory

### Lagrangian–linear strain relation

The Green-Lagrange strain tensor η and the linear strain ε are related by
(Landau & Lifshitz, §1):

```
η = ε + ½ ε²
```

This is solved iteratively to recover ε from a specified η:

```
ε_{n+1} = η − ½ εₙ²     (converges in < 20 steps for η ≲ 0.1)
```

The deformation matrix is `D = I + ε` and the deformed lattice vectors are:

```
R' = (D · Rᵀ)ᵀ
```

Strain values should stay in the harmonic regime (η ≲ 5–8 %).  Larger strains
can trigger phase transitions and invalidate the quadratic E(η) fit.

---

### Voigt shear convention

The Voigt shear component is twice the tensor component:

```
e₄ = 2η₂₃,   e₅ = 2η₁₃,   e₆ = 2η₁₂
```

So to apply a tensor shear of amplitude η, the Voigt vector carries **2η** in
the shear slot.  This is the standard convention used by ElaStic
(Golesorkhtabar & Pavone, CPC 2013) and the exciting reference implementation.

**Post-processing:** for a pure xy-shear run (code 6):

```
tensor η₁₂ = η
E(η) = E₀ + 2 V₀ C₄₄ η²
⇒  C₄₄ = d²E/dη² / (4 V₀)
```

---

### Volume conservation (cubic codes 1 and 2)

The [2,2] element of the strain tensor is set to:

```
η₃₃ = 1/(1 − η²) − 1
```

This enforces exact volume conservation for all η.  Proof for code 1 with
diagonal (η, −η, δ):

```
det(I + η) = (1 + η)(1 − η)(1 + δ) = (1 − η²)(1 + δ) = 1
⇒  δ = 1/(1 − η²) − 1
```

---

### Cubic elastic constants

Each deformation code isolates one Cᵢⱼ combination directly from the curvature
of E(η), without solving a linear system.

**Code 0** — volumetric (η, η, η, 0, 0, 0):
```
E(η) = E₀ + (9/2) V₀ B₀ η²   ⇒   B₀ = d²E/dη² / (9 V₀)
```

**Code 1** — orthorhombic, volume-conserving (η, −η, 1/(1−η²)−1, 0, 0, 0):
```
E(η) = E₀ + V₀ (C₁₁ − C₁₂) η²   ⇒   C₁₁ − C₁₂ = d²E/dη² / (2 V₀)
```

**Code 2** — monoclinic shear, volume-conserving (0, 0, 1/(1−η²)−1, 0, 0, 2η):
```
E(η) = E₀ + 2 V₀ C₄₄ η²   ⇒   C₄₄ = d²E/dη² / (4 V₀)
```

---

### Stress tensor

VASP stores stress in `vasprun.xml` in kbar (compression-positive).  `ase`
returns it in eV/Å³ with the tension-positive convention.

Unit conversion: **1 eV/Å³ = 160.21766208 GPa**.

Hydrostatic pressure and deviatoric stress:

```
P = −Tr(σ) / 3          (positive P = compression)
s = σ + P·I             (deviatoric; Tr(s) = 0)
```

---

## References

[1] exciting code — energy vs strain:
    http://exciting-code.org/nitrogen-energy-vs-strain-calculations

[2] ElaStic — Golesorkhtabar, Pavone et al., *Comput. Phys. Commun.* **184**,
    1861 (2013). https://doi.org/10.1016/j.cpc.2013.03.010

[3] Materials Project — elasticity calculations:
    https://wiki.materialsproject.org/Elasticity_calculations

[4] ELASTIC code (stress–strain approach):
    https://elastic.readthedocs.io/en/stable/

[5] IRelast — https://doi.org/10.1016/j.jallcom.2017.10.139

[6] L.D. Landau, E.M. Lifshitz, *Theory of Elasticity*, Elsevier (1986).
    ISBN: 0750626330
