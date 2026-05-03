# VASP-MDML-DataExtraction

This folder contains utility scripts to extract, convert and analyze VASP and ML-AIMD/MD datasets, prepare input/output for machine-learned interatomic potentials (MLIPs), and generate diagnostics such as potential energy surface (PES) maps and MSD/diffusion plots.

## Overview of scripts

- `1_MPtrj_Dataset_Query_and_Parsing_Criteria.py`
  - Query Materials Project tasks and main entries using mp-api and perform parsing/quality checks for MD/trajectory frames.
  - Implements energy convergence checks, GGA/GGA+U matching, and a uniqueness filter (structure matching with energy thresholds) to select unique trajectory frames.
  - Inputs: MP API key (set `MP_API_KEY`), material search criteria in the script.
  - Outputs: Prints matching/validation status for frames; can be adapted to save selected frames.

- `vaspouth5_to_extxyz.py`
  - Convert `vaspout.h5` (py4vasp) outputs into an Extended XYZ (`.extxyz`) trajectory for use with ASE and MD-analysis tools.
  - Writes lattice, per-atom positions, forces and stress and uses batching for large trajectories.
  - Inputs: `vaspout.h5` (py4vasp Calculation.from_file) and constants defined in the script.
  - Outputs: `EXTXYZ_NAME` (set in script) containing the trajectory in extended xyz format.

- `get_MLAIMDVASP_to_EXTXYZ.py`
  - Similar to the above but targeted at ML-AIMD data using py4vasp API interface; extracts energy/forces/stress and writes to an `.extxyz` file (`mlmd_from_outcar.extxyz`).
  - Inputs: `vaspout.h5` produced by py4vasp-compatible runs.
  - Outputs: `mlmd_from_outcar.extxyz` (buffered write)

- `generate_PESPlot_from_DFTData.py` and `get_PESPlot_DFTData.py`
  - Read `.extxyz` structures, compute pair-distance histograms as simple structural descriptors and reduce dimensionality with UMAP.
  - Interpolate energies on a 2D grid and output: smooth 3D PES (`smooth_PES_umap.png`), contour map (`contour_lines_umap.png`) and an interactive HTML (`smooth_PES_interactive.html`).
  - Inputs: Collections of `.extxyz` files (the script expects paths configured via `base_dirs` and file patterns such as `CONTCAR_mace.extxyz`).
  - Outputs: `train_data.extxyz` (intermediate) and the PES plots.

- `get_PureElementEnergyRef.py`
  - Compute per-element energy convergence vs vacuum size using `vasprun.xml` last-step energies.
  - Compares calculated energies to a MACE-MP reference energy table and produces tables and plots:
    - `energies.png`, `energies_convergence_subplots.png`, `energies_delta_log.png`.
  - Inputs: directory structure with vacuum subfolders and `vasprun.xml` files.
  - Outputs: CSV/console table and PNG plots.

- `plot_MSD_from_VASPKIT.py`
  - Read `MSD.dat` (and optionally `DIFFUSION_COEFFICIENT.dat`) produced by VASPKIT and produce interactive Plotly HTML plots (`msd_plot.html`, `diffusion_plot.html`).
  - Inputs: `MSD.dat`, `DIFFUSION_COEFFICIENT.dat`.

- `get_plot_MSD_GPU.py`
  - GPU-accelerated MSD computation using CuPy to parse `XDATCAR` and compute MSD on the GPU, then save `msd_data.csv` and `msd_plot.html`.
  - Inputs: `XDATCAR` (VASP trajectory), requires CuPy installed and a CUDA-capable GPU.

- `get_struct_from_Alexandria.py`
  - Downloader for the Alexandria dataset (PBE) hosted at alexandria.icams.rub.de; scrapes `*.bz2` files and downloads them in parallel using ThreadPoolExecutor.
  - Outputs: downloads into `complete_database/`.

- `test_model_from_SevenNet.py`
  - Evaluation utilities for the SevenNet potential (SevenNetCalculator). Evaluates parity (DFT vs base vs fine-tuned model) for energy/forces/stress, generates parity plots and EOS comparisons.
  - Inputs: `data/evaluation` directory (expects `test_md.extxyz`, `eos.extxyz`, and model checkpoints like `checkpoint_fine_tuned.pth`).
  - Outputs: parity plots (PNG), EOS `.extxyz` outputs under `eos/`.


## Dependencies

Most scripts use common scientific Python libraries. Typical dependencies include:

- Python 3.8+
- ase
- numpy
- scipy
- scikit-learn
- umap-learn
- joblib
- matplotlib
- pandas
- plotly
- tqdm
- py4vasp (for reading `vaspout.h5`)
- monty (optional for some utilities)
- py-matgen / mp-api (for Materials Project scripts)
- cupy (for GPU MSD script)
- sevenn (for SevenNet scripts)
- requests, beautifulsoup4 (for Alexandria downloader)

Install with pip or conda, e.g.:

pip install ase numpy scipy scikit-learn umap-learn joblib matplotlib pandas plotly tqdm py4vasp monty pymatgen mp-api cupy-baremetal beautifulsoup4 requests

(Adjust for conda and GPU-specific CuPy packages.)


## Usage notes and examples

- Convert py4vasp `vaspout.h5` to extxyz (batch-friendly):
  - Ensure `py4vasp` can open `vaspout.h5` in current directory, then run:
    python vaspouth5_to_extxyz.py

- Create `.extxyz` from ML-AIMD/outcar-style h5 files:
    python get_MLAIMDVASP_to_EXTXYZ.py

- Generate PES maps from a folder of `.extxyz` or `CONTCAR_mace.extxyz` files:
    - Edit `base_dirs` in `generate_PESPlot_from_DFTData.py` or `get_PESPlot_DFTData.py` to point to your data and run.

- GPU-accelerated MSD from XDATCAR:
    python get_plot_MSD_GPU.py

- Materials Project trajectory filtering (requires MP API key):
    - Set `MP_API_KEY` in `1_MPtrj_Dataset_Query_and_Parsing_Criteria.py` and adapt search criteria, then run.


## Output conventions

- Extended XYZ files include lattice, per-atom forces and stress and are compatible with ASE and many ML/MD analysis tools.
- Plots are saved as high-resolution PNGs and interactive HTMLs where applicable.


## Suggestions / Next steps

- Add small example data or a `data/` pointer for quick local testing of plotting scripts.
- Add command-line argument parsing (argparse) to each script to configure inputs/outputs without editing the source.
- Add unit tests / smoke tests (e.g., tiny example trajectory) to validate conversions.


---
