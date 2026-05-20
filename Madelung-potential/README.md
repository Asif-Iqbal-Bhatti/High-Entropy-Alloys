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

### 3. Units and Conversions
To provide results in standard units (**Volts** and **electron-Volts**), the script uses the following conversion:
The conversion factor from internal units ($e/\mathring{A}$) to Volts is:

$$ 1 \text{ unit} = \frac{e}{4\pi\epsilon_0 \cdot 1 \mathring{A}} \approx 14.3996 \text{ V} $$

---

## Usage
1. Ensure your structure file is named `CONTCAR` or `POSCAR`.
2. Append your partial charges as a **4th column** in the coordinate section.
