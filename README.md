# High-Entropy-Alloys — Computational Toolkit

A collection of Python (and some C++/Fortran) scripts for first-principles and machine-learning-based analysis of High-Entropy Alloys (HEAs). The tools cover the full workflow from structure generation and thermodynamics through defect mechanics, machine-learned potentials, and deep-learning property prediction.

Scripts accompanying the peer-reviewed publication:
> [doi.org/10.1016/j.commatsci.2024.113196](https://doi.org/10.1016/j.commatsci.2024.113196)

---

## Repository layout

```
High-Entropy-Alloys/
├── src/                              # Core structure and thermodynamics scripts
├── MLIP_MonteCarlo/                  # Monte Carlo with machine-learned potentials
├── VASP_MLFF/                        # VASP machine-learned force-field utilities
├── VASP-MDML-DataExtraction/         # VASP → Extended-XYZ / ML-dataset preparation
├── MaterialsPropertyAnalysis_ML_DNN/ # ML / deep-learning property prediction
├── Elastic_tensor_from_Lammps/       # Elastic tensor extraction via LAMMPS
├── surface_states_alloys/            # Surface energy calculations
└── utilities/                        # Plotting, formation energy, SRO, dislocations, …
```

---

## Modules

### `src/` — Structure generation & core thermodynamics

| Script | Purpose |
|--------|---------|
| `SQS_structure_gen.py` | Generates Special Quasi-random Structures (SQS) for multi-component HEAs |
| `mcsqs2poscar.py` / `sqs2poscar_modified_AIB.cpp` | Converts MCSQS output to VASP POSCAR format |
| `get_mcsqs_lattices.py` | Prepares BCC/FCC/HCP lattice inputs for MCSQS |
| `PMSROCalculatorPyScale2.py` | Computes PM-SRO and WC-SRO short-range order parameters across neighbour shells |
| `get_ConfigurationEntropy_HEA.py` | Calculates configurational entropy (ΔS_conf = −k_B Σ c_i ln c_i) from CIF files, with temperature sweep and CSV export |
| `GSFE_BCC.py` / `GSFE.sh` | Generalized Stacking Fault Energy (GSFE) calculations for BCC crystals |
| `BCCScrew_dislocation.py` / `BCCScrew_dislocation.f90` | Screw dislocation core analysis for BCC lattices |
| `lattice_mismatch_measure.py` | Quantifies lattice mismatch between alloy components |
| `scalingPOSCAR.py` / `supercell_rand.py` | POSCAR scaling and random supercell construction |
| `set_orderTypes_POSCAR.py` / `set_shake_atom_lattice.py` | Atomic ordering and thermal displacement tools |
| `read_POSCAR.cpp` | Lightweight C++ POSCAR reader |

---

### `MLIP_MonteCarlo/` — Monte Carlo with machine-learned potentials

`MLIP_MonteCarlo_HEA.py` — Metropolis Monte Carlo driven by a GRACE machine-learned interatomic potential. Performs random atomic substitutions on a supercell, accepts/rejects moves based on energy at a target temperature, and saves the lowest-energy configuration found.

---

### `VASP_MLFF/` — VASP machine-learned force-field utilities

| Script | Purpose |
|--------|---------|
| `get_data.py` | Extracts training frames from VASP ML-FF runs |
| `get_MLData_py4vasp.py` | Reads VASP output via py4vasp for ML dataset construction |
| `get_VASP_ML_ERR_Analysis.py` | Analyses force/energy errors of the VASP ML-FF against DFT reference |

---

### `VASP-MDML-DataExtraction/` — VASP → ML-ready datasets

| Script | Purpose |
|--------|---------|
| `vaspouth5_to_extxyz.py` | Converts `vaspout.h5` trajectories to Extended XYZ |
| `get_MLAIMDVASP_to_EXTXYZ.py` | Extracts AIMD frames from VASP output to Extended XYZ |
| `1_MPtrj_Dataset_Query_and_Parsing_Criteria.py` | Queries and filters the MPtrj dataset (energy convergence, GGA/GGA+U matching) |
| `get_PureElementEnergyRef.py` | Retrieves DFT reference energies for pure elements |
| `generate_PESPlot_from_DFTData.py` / `get_PESPlot_DFTData.py` | Potential energy surface visualisation via UMAP dimensionality reduction |
| `get_plot_MSD_GPU.py` / `plot_MSD_from_VASPKIT.py` | GPU-accelerated (CuPy) and standard MSD analysis |
| `get_struct_from_Alexandria.py` | Downloads structures from the Alexandria database |
| `test_model_from_SevenNet.py` | Evaluates a SevenNet MLIP against DFT data |

---

### `MaterialsPropertyAnalysis_ML_DNN/` — Machine learning & deep learning

| Script | Purpose |
|--------|---------|
| `ML_DNN_with_TensorFlowKeras/DeepNeuralNetwrok_TF_HEA.py` | Deep neural network (TF/Keras) for HEA property regression |
| `ML_DNN_with_TensorFlowKeras/LinearRegression_TF_HEA.py` | TensorFlow linear regression baseline |
| `ML_DNN_with_TensorFlowKeras/get_gaussianProcessRegression.py` | Gaussian process regression for property prediction |
| `pred_short_range_order_pytorch.py` | PyTorch neural network trained on Materials Project data to predict SRO / alloy properties |

---

### `Elastic_tensor_from_Lammps/`

`elastic_tensor_lammps.py` — Automates extraction and symmetrisation of the full elastic stiffness tensor from LAMMPS deformation runs.

---

### `surface_states_alloys/`

| Script | Purpose |
|--------|---------|
| `1_generate_surface_pmg.py` | Generates surface slabs using Pymatgen |
| `get_surface_energy_using_pymatgen.py` | Computes surface energy (J/m²) from VASP CONTCAR/OUTCAR files relative to a bulk reference |

---

### `utilities/` — Supporting tools

| Sub-directory / script | Purpose |
|------------------------|---------|
| `HEA_composition.py` | Composition parsing and weight-fraction ↔ atom-fraction conversion |
| `configuration_entropy.py` | Standalone configurational entropy calculator |
| `coordination_analysis.py` | Coordination number analysis from trajectory data |
| `Poscar_colorPlotter.py` | Colour-coded POSCAR structure visualisation |
| `get_InitialRelax_from_ASE.py` | Initial structure relaxation via ASE |
| `get_dispersion.py` | Phonon/electronic dispersion utilities |
| `ClusterExpansion_tools/icet_ClusterExpansion.py` | Cluster expansion with ICET |
| `ConversionData_tools/` | Unit/format converters (DAT → JSON, wt% → at%) |
| `DislocationDipole_tools/` | Dislocation dipole elastic energy, Frenkel-Kontorova model, core energy |
| `FormEnergy_tools/calculate_formation_energy.py` | Formation energy calculation from DFT total energies |
| `MCSQS_tools/` | Parallel SQS generation, SRO via pyscal, MCSQS file numbering |
| `MonteCarlo_tools/Metropolis_Monte_carlo.py` | Metropolis Monte Carlo reference implementation |
| `PhaseDiagram_tools/phase_diagram_for_HEA.py` | HEA phase diagram construction |
| `plot_tools_HEA/` | Vitek maps, box plots, bar charts, spider plots, 3-D spline fitting, PDOS plotting, Fresnel rendering |

---

## Dependencies

| Category | Packages |
|----------|----------|
| Core scientific | `numpy`, `scipy`, `matplotlib`, `pandas` |
| Atomistic | `ase`, `pymatgen`, `pyscal` |
| VASP interface | `py4vasp` |
| Machine learning | `scikit-learn`, `tensorflow` / `keras`, `torch` |
| MLIPs | `grace` (GRACE potential), `sevenn` (SevenNet) |
| Dimensionality reduction | `umap-learn`, `plotly` |
| GPU acceleration | `cupy` (optional, for MSD on CUDA hardware) |
| Other | `icet` (cluster expansion), `atomsk` (structure generation, external binary) |

Python 3.8+ is required throughout.

---

## Citation

If you use these scripts in your research, please cite:

> Bhatti, A. I. *et al.*, *Computational Materials Science*, 2024.
> [https://doi.org/10.1016/j.commatsci.2024.113196](https://doi.org/10.1016/j.commatsci.2024.113196)

---

## License

Distributed under the **GNU General Public License v3.0**. See [`LICENCE`](LICENCE) for details.

> **Disclaimer:** Scripts are provided as-is with no warranty. Please review each script carefully before use in a research workflow.
