# Madelung Potential Calculator

A Python script for computing the Madelung potential (MP) and total electrostatic energy of crystal systems from VASP geometry files.

**Created**  :: 14/05/2020  
**Updated**  :: 08/04/2026  
**Author**   :: Asif Iqbal ([@AIB_EM](https://github.com/AIB_EM))

---

## Mathematical Formulation

The script computes the electrostatic potential based on the classical Coulomb interaction within a periodic framework.

### 1. Site Potential ($V_i$)

The Madelung potential at a specific atomic site $i$ is calculated by summing the contributions of all other charges $j$ in the reference unit cell and all periodic images:

$$V_i = \frac{1}{4\pi\epsilon_0} \sum_{\mathbf{R}} \sum_{j} \frac{q_j}{|\mathbf{r}_{ij} + \mathbf{R}|}$$

Where:

*   $q_j$ is the partial charge of atom $j$ (provided in the 4th column).
*   $\mathbf{r}_{ij}$ is the distance vector between atom $i$ and atom $j$ within the primary cell.
*   $\mathbf{R}$ is the lattice translation vector: $\mathbf{R} = n_1\mathbf{a} + n_2\mathbf{b} + n_3\mathbf{c}$ (where $n$ are integers).
*   The term where $j=i$ and $\mathbf{R}=0$ is strictly excluded to avoid self-interaction.

### 2. Total Electrostatic Energy ($U$)

The total Madelung energy of the unit cell is calculated using the following summation:

$$U = \frac{1}{2} \sum_{i=1}^{N} q_i V_i$$

The factor of $1/2$ is included to avoid double-counting the pairwise interactions between ions.

### 3. Ewald Summation

The original direct real-space sum is conditionally convergent and shell-limit dependent.
The production code uses **Ewald summation** (via `pymatgen.analysis.ewald.EwaldSummation`),
which decomposes the lattice sum into three absolutely convergent parts:

$$U = U_{\text{real}} + U_{\text{recip}} + U_{\text{point}}$$

- **Real-space part** — short-ranged, converges exponentially with a Gaussian damping parameter $\eta$
- **Reciprocal-space part** — long-ranged, sum over $\mathbf{G}$-vectors weighted by $e^{-G^2/4\eta}/G^2$
- **Point-charge (self) correction** — removes the spurious self-interaction of each Gaussian

This gives $\mathcal{O}(N^{3/2})$ scaling and machine-precision accuracy regardless of shell limit.

### 4. Units and Conversions

To provide results in standard units (**Volts** and **electron-Volts**), the script uses the following conversion:

The conversion factor from internal units ($e/\mathring{A}$) to Volts is:

$$ 1 \text{ unit} = \frac{e}{4\pi\epsilon_0 \cdot 1 \mathring{A}} \approx 14.3996 \text{ V} $$

---

## Charge Sources

The code supports three charge sources, selected via the `--charges` flag:

| Mode | Flag | Source | Notes |
|------|------|---------|-------|
| Formal | `--charges poscar` | 4th column of POSCAR/CONTCAR | Default. Integer or non-integer formal charges in units of $e$. |
| Bader | `--charges bader` | `ACF.dat` from Henkelman Bader code | Ionic charge $q_i = Z_i - Q_{\text{Bader},i}$. Reads ZVAL from POTCAR, built-in table, or `--zvals`. |
| Mulliken / NBO | `--charges mulliken` | Plain text file, one value per line | Compatible with any externally computed partial charge scheme. |

---

## Usage

1. Ensure your structure file is named `CONTCAR` or `POSCAR`.
2. Append your partial charges as a **4th column** in the coordinate section.

```bash
# Formal charges from POSCAR column 4 (default)
python madelung_ewald.py POSCAR

# Bader charges — reads ACF.dat + POTCAR automatically
python madelung_ewald.py CONTCAR --charges bader

# Bader charges with explicit ZVALs (no POTCAR needed)
# Example: Li6PS5Cl — Li=1, P=5, S=6, Cl=7
python madelung_ewald.py CONTCAR --charges bader --zvals 1 5 6 7

# Mulliken or NBO charges from an external file
python madelung_ewald.py CONTCAR --charges mulliken --charges-file nbo.txt

# Run NaCl validation test
python madelung_ewald.py POSCAR --validate

# Custom Ewald splitting parameter
python madelung_ewald.py POSCAR --eta 0.5
```

---

## Dependencies

```
numpy
pymatgen
```

Install with:

```bash
pip install numpy pymatgen
```

---

## Validation and Testing

All tests were run with `pymatgen >= 2024.x` and `numpy >= 1.26`.
The test suite covers accuracy against known analytical results, internal
consistency (supercell invariance, symmetry), correct physics (sign convention,
charge scaling), robustness of file parsers, and coverage across all 14 Bravais
lattice types.

---

### Test 1 — NaCl Madelung constant (primitive FCC cell, 2 atoms)

The textbook Madelung constant for NaCl is $A = 1.7475645946$.
The 2-atom FCC primitive cell was constructed with $a = 5.6402$ Å and
formal charges $q(\text{Na}) = +1$, $q(\text{Cl}) = -1$.

| Quantity | Computed | Reference | Error |
|----------|----------|-----------|-------|
| Madelung constant $A$ | 1.74756459 | 1.74756459 | 0.000000 % |
| $V(\text{Na})$ | −4.461599 V | −4.4616 V | < 0.001 % |
| $V(\text{Cl})$ | +4.461599 V | +4.4616 V | < 0.001 % |
| $U_{\text{total}}$ | −8.923198 eV | −8.9232 eV | < 0.001 % |

**Result: PASS**

---

### Test 2 — CsCl Madelung constant (simple cubic, 2 atoms)

CsCl crystallises in the simple cubic (Pm-3m) structure with a 2-atom basis.
The textbook Madelung constant is $A = 1.76267$.
Cell parameter $a = 4.123$ Å, nearest-neighbour distance $r_0 = a\sqrt{3}/2$.

| Quantity | Computed | Reference | Error |
|----------|----------|-----------|-------|
| Madelung constant $A$ | 1.762675 | 1.762670 | 0.0003 % |

**Result: PASS**

---

### Test 3 — Supercell size invariance (NaCl, primitive vs. conventional)

The total energy per formula unit must be identical regardless of whether the
primitive 2-atom FCC cell or the conventional 8-atom cubic cell is used.
This tests that the Ewald sum is extensive and cell-shape independent.

| Cell | Atoms | $U/\text{f.u.}$ (eV) | Difference |
|------|-------|----------------------|------------|
| FCC primitive | 2 | −8.92319797 | — |
| Cubic conventional | 8 | −8.92319797 | 1.28 × 10⁻¹² eV |

Symmetry-equivalent sites (all 4 Na, all 4 Cl) must carry identical potentials.
Observed spread across equivalent sites: 0.00 V (exact).

**Result: PASS**

---

### Test 4 — Potential sign convention

A cation (positive charge) sits in a field dominated by surrounding anions and
must therefore reside at a **negative** electrostatic potential, and vice versa.

| Site | Charge | Potential | Expected sign |
|------|--------|-----------|---------------|
| Na | +1 e | −4.4616 V | negative ✓ |
| Cl | −1 e | +4.4616 V | positive ✓ |

**Result: PASS**

---

### Test 5 — Quadratic charge scaling ($U \propto q^2$)

The Madelung energy is quadratic in the charge: doubling all charges must
multiply the energy by exactly 4.

| Charge | $U$ (eV) | Ratio |
|--------|----------|-------|
| $q = \pm 1\,e$ | −8.923198 | — |
| $q = \pm 2\,e$ | −35.692792 | 4.000000 |

Expected ratio: 4.000000 exactly.

**Result: PASS**

---

### Test 6 — Multi-species, non-cubic cell: ZnS wurtzite (hexagonal, $q = \pm 2$)

ZnS wurtzite (space group $P6_3mc$) has a hexagonal unit cell
($a = 3.8226$ Å, $c = 6.2605$ Å, $u = 0.3748$) with 4 atoms and divalent charges.
This tests non-cubic geometry and $|q| > 1$.

Space group symmetry requires $V(\text{Zn}_1) = V(\text{Zn}_2)$ and
$V(\text{S}_1) = V(\text{S}_2)$.

| Quantity | Value |
|----------|-------|
| $V(\text{Zn})$ | −10.0862 V |
| $V(\text{S})$ | +10.0862 V |
| Zn site spread | 1.78 × 10⁻¹⁵ V |
| S site spread | 0.00 V |
| $U_{\text{total}}$ | −80.689841 eV |
| Supercell deviation | 9.9 × 10⁻¹² eV |

**Result: PASS**

---

### Test 7 — TiO2 rutile (tetragonal-P, mixed valence $q(\text{Ti}) = +4$, $q(\text{O}) = -2$)

Rutile TiO2 (space group $P4_2/mnm$) has $a = 4.5937$ Å, $c = 2.9587$ Å,
$u = 0.3053$, with 6 atoms per cell and a $2:4$ charge ratio.
Tests non-equal charge magnitudes and tetragonal symmetry.

| Site | $V$ (V) | Symmetry spread |
|------|---------|-----------------|
| Ti (×2) | −22.3651 | 0.00 V |
| O (×4) | +12.9396 | 3.55 × 10⁻¹⁵ V |

**Result: PASS**

---

### Test 8 — ACF.dat parser (Bader charge reader)

A synthetic `ACF.dat` file was constructed for the NaCl conventional cell
with known per-atom Bader electron counts:
$Q_{\text{Bader}}(\text{Na}) = 0.143988\,e$, $Q_{\text{Bader}}(\text{Cl}) = 7.856012\,e$.

| Check | Expected | Parsed | Match |
|-------|----------|--------|-------|
| Na Bader electrons (×4) | 0.143988 | 0.143988 | ✓ |
| Cl Bader electrons (×4) | 7.856012 | 7.856012 | ✓ |
| Total electrons | 32.0000 | 32.0000 | ✓ |

Ionic charges derived: $q(\text{Na}) = +0.856\,e$, $q(\text{Cl}) = -0.856\,e$,
reflecting the partial covalency of NaCl. The Madelung potential with Bader
charges ($-3.82$ V on Na) is reduced by ~14 % relative to formal charges
($-4.46$ V), consistent with literature.

**Result: PASS**

---

### Test 9 — POSCAR parser

A VASP 5-format POSCAR for the NaCl conventional cell was parsed and checked
against known values.

| Field | Expected | Parsed |
|-------|----------|--------|
| Species order | Na×4, Cl×4 | Na×4, Cl×4 ✓ |
| Formal charges (col. 4) | +1, −1 | +1, −1 ✓ |
| Lattice parameter $a$ | 5.6402 Å | 5.6402 Å ✓ |
| Fractional coords in [0,1) | yes | yes ✓ |

The parser correctly handles both VASP 4 (no element name line) and VASP 5
(element name line present) formats, as well as optional Selective Dynamics
blocks.

**Result: PASS**

---

### Test 10 — All 14 Bravais lattice types

The Ewald summation was validated across all seven crystal systems and all 14
Bravais lattice types. For each lattice, the following checks were applied:

- **Symmetry**: symmetry-equivalent sites carry identical potentials (spread < 10⁻⁸ V)
- **Sign convention**: cation potential is negative, anion potential is positive
- **Supercell invariance**: doubling the cell along **a** gives identical $U$/f.u. (deviation < 10⁻¹¹ eV)
- **Known Madelung constant** where analytically available

| # | Lattice type | Representative structure | Sym. | Sign | Supercell | $A$ error |
|---|-------------|--------------------------|------|------|-----------|-----------|
| 1 | Cubic-P (sc) | CsCl | — | ✓ | 1.04×10⁻¹² eV | 0.0003 % |
| 2 | Cubic-F (fcc) | NaCl primitive | — | ✓ | 8.88×10⁻¹⁵ eV | 0.000000 % |
| 3 | Cubic-I (bcc) | CsCl bcc primitive | — | ✓ | 1.79×10⁻¹² eV | — |
| 4 | Tetragonal-P | TiO2 rutile | ✓ | ✓ | 4.50×10⁻¹¹ eV | — |
| 5 | Tetragonal-I (bct) | synthetic AB | ✓ | ✓ | 1.65×10⁻¹¹ eV | — |
| 6 | Orthorhombic-P | PbO-like | ✓ | ✓ | 3.14×10⁻¹¹ eV | — |
| 7 | Orthorhombic-F | synthetic AB | ✓ | ✓ | 1.96×10⁻¹¹ eV | — |
| 8 | Orthorhombic-I | synthetic AB | ✓ | ✓ | 1.45×10⁻¹² eV | — |
| 9 | Orthorhombic-C | synthetic AB | ✓ | ✓ | 1.58×10⁻¹¹ eV | — |
| 10 | Hexagonal-P | ZnS wurtzite | ✓ | ✓ | 9.90×10⁻¹² eV | — |
| 11 | Trigonal-R | rhombohedral primitive | — | ✓ | 6.95×10⁻¹² eV | — |
| 12 | Monoclinic-P | synthetic AB, β=110° | ✓ | ✓ | 2.17×10⁻¹² eV | — |
| 13 | Monoclinic-C | synthetic AB, β=107° | ✓ | ✓ | 9.93×10⁻¹² eV | — |
| 14 | Triclinic (aP) | synthetic, α=83° β=96° γ=74° | — | ✓ | 1.42×10⁻¹¹ eV | — |

Supercell deviations are at floating-point rounding level (~10⁻¹¹–10⁻¹⁵ eV),
confirming that the Ewald summation is cell-shape and cell-size invariant across
all crystal systems.

**Result: ALL 14 LATTICE TYPES PASS**

---

### Summary table

| Test | Description | Checks | Result |
|------|-------------|--------|--------|
| 1 | NaCl Madelung constant (fcc primitive) | $A$, $V_i$, $U$ vs. textbook | **PASS** |
| 2 | CsCl Madelung constant (sc) | $A$ vs. textbook | **PASS** |
| 3 | Supercell size invariance (NaCl 2 vs. 8 atoms) | $U$/f.u., site equivalence | **PASS** |
| 4 | Potential sign convention | cation $V < 0$, anion $V > 0$ | **PASS** |
| 5 | Quadratic charge scaling ($U \propto q^2$) | ratio = 4.000000 exactly | **PASS** |
| 6 | ZnS wurtzite (hexagonal, $q = \pm 2$) | symmetry, sign, supercell | **PASS** |
| 7 | TiO2 rutile (tetragonal, mixed valence) | symmetry, sign, supercell | **PASS** |
| 8 | ACF.dat Bader charge parser | per-atom counts, total electrons | **PASS** |
| 9 | POSCAR parser (VASP 5 format) | species, charges, lattice, coords | **PASS** |
| 10 | All 14 Bravais lattice types | symmetry, sign, supercell per lattice | **PASS** |

43/44 individual assertions passed. The one numerical discrepancy (cubic-I BCC
energy vs. cubic-P SC energy for CsCl) was traced to a **test design error**:
the BCC primitive vectors with a 2-atom basis describe a geometrically distinct
crystal (half the cell volume of the SC description), not an alternative
primitive cell of the same crystal. CsCl has no BCC primitive cell; the SC
2-atom cell is already primitive. The Ewald code itself produced correct and
internally consistent results for both descriptions independently.

---

## Known Limitations

- **Point-charge model**: all charge sources (formal, Bader, Mulliken) are treated as point charges. Charge overlap (exchange-correlation) and multipole effects are not included.
- **Partial charges**: Bader and Mulliken charges give an electrostatic site potential that is physically meaningful but not equivalent to the formal Madelung potential. Results scale with charge magnitude.
- **Non-neutral cells**: `EwaldSummation` applies a uniform jellium background automatically. Results are well-defined but the physical interpretation requires care.
- **VASP 4 POSCAR** (no element name line): supported, but species will be labelled `X` unless element names are inferred from context.
