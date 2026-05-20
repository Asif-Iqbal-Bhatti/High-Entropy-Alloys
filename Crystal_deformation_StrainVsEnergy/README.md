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

Voigt convention: normal components eᵢ = η; shear components eᵢ = **2η** so
that the tensor off-diagonal = eᵢ/2 = η.  This is consistent with ElaStic and
exciting (see [Voigt shear convention](#voigt-shear-convention) below).

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

Codes 1 and 2 enforce exact volume conservation for all η (see
[Volume conservation](#volume-conservation) below).

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

See [Stress sign convention](#stress-sign-convention) for details on the VASP
kbar → ase eV/Å³ conversion and the sign flip.

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

All 80 tests are self-contained (no VASP/exciting files required):

```
tests/test_strain_core.py          core numerical kernels
tests/test_lagrangian_strain.py    POSCAR I/O, deformation codes, shear convention
tests/test_cubic_strain.py         cubic η matrices, volume conservation
tests/test_stress_pressure.py      pressure and deviatoric stress
```

---

## Theory and implementation notes

### Lagrangian–linear strain relation

The Green-Lagrange strain tensor η and the linear (infinitesimal) strain ε are
related by (Landau & Lifshitz, §1):

```
η = ε + ½ ε²
```

This is solved iteratively to recover ε from a specified η:

```
ε_{n+1} = η − ½ εₙ²     (converges in < 20 steps for η ≲ 0.1)
```

The deformation matrix is then `D = I + ε`, and deformed lattice vectors are:

```
R' = (D · Rᵀ)ᵀ
```

This relation is confirmed exact by Landau & Lifshitz and verified numerically:
the solver residual `‖η − (ε + ½ε²)‖` reaches machine precision (< 10⁻¹⁶)
and the result matches the closed form `ε_diag = √(1 + 2η) − 1` for diagonal η.

Strain values should stay in the harmonic regime (η ≲ 5–8 %).  Larger strains
can trigger phase transitions and invalidate the quadratic E(η) fit.

---

### Voigt shear convention

This is the most common source of confusion and factor-of-2 errors when
comparing codes.

**The Voigt shear component eₖ is twice the corresponding tensor component:**

```
e₄ = 2η₂₃,   e₅ = 2η₁₃,   e₆ = 2η₁₂
```

So to apply a tensor shear strain of amplitude η (e.g. for C₄₄), the Voigt
vector must carry **2η** in the shear slot, not η.

Both `lagrangian_strain.py` and `cubic_strain.py` follow this convention,
consistent with the ElaStic code (Golesorkhtabar & Pavone, CPC 2013) and the
exciting reference implementation.

**Consequence for post-processing:**

For a pure xy-shear deformation (code 6, Voigt = (0,0,0,0,0,2η)):

```
tensor η₁₂ = η
E(η) = E₀ + 2 V₀ C₄₄ η²  +  O(η⁴)
⇒  C₄₄ = (1/V₀) · d²E/dη² / 4
```

Your analysis script must divide the fitted second-order coefficient by **4V₀**
(or equivalently divide the leading coefficient of the η² polynomial by 2V₀,
since 2 V₀ C₄₄ is the coefficient).  This is what ElaStic does.

> **Warning for future readers:** An earlier version of these scripts used
> `'E'` (Voigt = η) at shear positions instead of `'2'` (Voigt = 2η).
> This gave tensor shear = η/2, meaning C₄₄ was extracted directly without the
> factor of 4 — a different but internally consistent convention.  The scripts
> were updated to use the standard ElaStic/exciting convention so that the
> same post-processing workflow can be used unchanged.

---

### Volume conservation

For cubic codes 1 and 2, the [2,2] element of the Lagrangian strain tensor is
set to:

```
η₃₃ = 1/(1 − η²) − 1
```

This enforces **exact** volume conservation for all η (not just to second
order).  Proof for code 1 with diagonal (η, −η, δ):

```
det(I + η) = (1 + η)(1 − η)(1 + δ) = (1 − η²)(1 + δ) = 1
⇒  δ = 1/(1 − η²) − 1   (exact for all η < 1)
```

Verified numerically: `det(I + η_mat) = 1.000000000000` for η = 0.01, 0.05,
0.10 for both codes 1 and 2.

---

### Cubic elastic constants — deformation matrix derivation

The three cubic deformation codes are chosen to isolate one elastic constant
combination per run, avoiding the need to solve a linear system.

**Code 0** — volumetric strain (η, η, η, 0, 0, 0):
```
E(η) = E₀ + (9/2) V₀ B₀ η²   ⇒   d²E/dη² = 9 V₀ B₀
```
where B₀ = (3C₁₁ + 6C₁₂)/9 is the bulk modulus.

**Code 1** — orthorhombic, volume-conserving (η, −η, 1/(1−η²)−1, 0, 0, 0):
```
E(η) = E₀ + V₀ (C₁₁ − C₁₂) η²   ⇒   d²E/dη² = 2 V₀ (C₁₁ − C₁₂)
```

**Code 2** — monoclinic shear, volume-conserving (0, 0, 1/(1−η²)−1, 0, 0, 2η):
```
tensor η₁₂ = η
E(η) = E₀ + 2 V₀ C₄₄ η²   ⇒   d²E/dη² = 4 V₀ C₄₄
```

These match the IRelast/WIEN2k formulas (Ref [5]).  ElaStic uses a different set
of three deformation codes and then solves a 3×3 linear system; both strategies
yield identical Cᵢⱼ values.

---

### Stress sign convention

VASP stores the stress in `vasprun.xml` in **kbar** with the compression-positive
sign (i.e. it stores −σ in the physics convention).

`ase.io.read` returns the stress in **eV/Å³** with the tension-positive
convention (i.e. it has already flipped the VASP sign).  The scripts use the ase
convention throughout.

Unit conversion: **1 eV/Å³ = 160.21766208 GPa**.

The hydrostatic pressure and deviatoric stress are computed as:

```
P = −Tr(σ) / 3          (positive P = compression)
s = σ + P·I             (deviatoric; traceless: Tr(s) = 0)
σ = s − P·I             (reconstruction)
```

> **Warning for future readers:** An earlier version of the script computed
> `s = σ − P·I`, which has the wrong sign and gives a non-zero deviatoric
> for isotropic stress.  The correct formula is `s = σ + P·I` (equivalently,
> `s = σ − (Tr(σ)/3)·I`).  This was cross-verified: for isotropic
> σ = σ₀ I, the corrected code gives s = 0 exactly.

---

## References

[1] exciting code — energy vs strain:
    http://exciting-code.org/nitrogen-energy-vs-strain-calculations

[2] ElaStic: A tool for calculating second-order elastic constants from first
    principles — Golesorkhtabar, Pavone et al., *Comput. Phys. Commun.* **184**,
    1861 (2013). https://doi.org/10.1016/j.cpc.2013.03.010

[3] Materials Project — elasticity calculations:
    https://wiki.materialsproject.org/Elasticity_calculations

[4] ELASTIC code (stress–strain approach):
    https://elastic.readthedocs.io/en/stable/

[5] IRelast package — Voigt deformation codes for cubic systems:
    https://doi.org/10.1016/j.jallcom.2017.10.139

[6] L.D. Landau, E.M. Lifshitz, *Theory of Elasticity*, Elsevier (1986).
    ISBN: 0750626330
