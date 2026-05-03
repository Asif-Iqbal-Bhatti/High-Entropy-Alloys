#!/usr/bin/env python
# ----------------------------- Imports -----------------------------
import os
import numpy as np
import torch
import matplotlib.pyplot as plt
from ase.io import read, write
from ase.constraints import UnitCellFilter
from ase.optimize import LBFGS
from ase.calculators.singlepoint import SinglePointCalculator
from scipy.stats import gaussian_kde
from tqdm import tqdm
from copy import deepcopy
from sevenn.calculator import SevenNetCalculator
import warnings

# ----------------------------- Configuration -----------------------------
warnings.filterwarnings("ignore")
WORKING_DIR = './'
DATA_PATH = 'data'
TO_KBAR = 1602.1766208

# ----------------------------- Utility Functions -----------------------------
def flatten_arrays(arrays):
    return np.concatenate([a.reshape(-1,) for a in arrays])

def load_trajectory(path):
    return read(path, index=':')

def density_colored_scatter_plot(dft_energy, nnp_energy, dft_force, nnp_force, dft_stress, nnp_stress, title=None):
    unit = {"energy": "eV/atom", "force": r"eV/$\\rm{\\AA}$", "stress": "kbar"}
    modes = ['energy', 'force', 'stress']
    plt.figure(figsize=(18/2.54, 6/2.54))
    for num, (x, y) in enumerate(zip([dft_energy, dft_force, dft_stress], [nnp_energy, nnp_force, nnp_stress])):
        mode = modes[num]
        idx = np.random.choice(len(y), 1000) if len(y) > 1000 else list(range(len(y)))
        xsam = [x[i] for i in idx]
        ysam = [y[i] for i in idx]
        xy = np.vstack([x, y])
        xysam = np.vstack([xsam, ysam])
        zsam = gaussian_kde(xysam)
        z = zsam.pdf(xy)
        idx = z.argsort()
        x = [x[i] for i in idx]
        y = [y[i] for i in idx]
        z = [z[i] for i in idx]
        ax = plt.subplot(int(f'13{num+1}'))
        plt.scatter(x, y, c=z, s=4, cmap='plasma')
        mini = min(min(x), min(y))
        maxi = max(max(x), max(y))
        ran = (maxi-mini) / 20
        plt.plot([mini-ran, maxi+ran], [mini-ran, maxi+ran], color='grey', linestyle='dashed')
        plt.xlim(mini-ran, maxi+ran)
        plt.ylim(mini-ran, maxi+ran)
        plt.xlabel(f'DFT {mode} ({unit[mode]})')
        plt.ylabel(f'MLP {mode} ({unit[mode]})')
        ax.set_aspect('equal')
        if title:
            ax.set_title(f'{title} {mode}')
    plt.tight_layout()
    plt.savefig(f"parity_plot_DFT_ML-{title}.png", bbox_inches="tight", dpi=300)
    plt.show()

# ----------------------------- Parity Evaluation -----------------------------
def evaluate_parity(traj, calc_fine, calc_base):
    dft_energy, dft_forces, dft_stress = [], [], []
    ft_energy, ft_forces, ft_stress = [], [], []
    base_energy, base_forces, base_stress = [], [], []

    for atoms in tqdm(traj, leave=False, desc='Getting data'):
        atoms.calc = calc_fine
        ft_energy.append(atoms.get_potential_energy() / len(atoms))
        ft_forces.append(atoms.get_forces())
        ft_stress.append(-1 * atoms.get_stress(voigt=False) * TO_KBAR)

        atoms.calc = calc_base
        base_energy.append(atoms.get_potential_energy() / len(atoms))
        base_forces.append(atoms.get_forces())
        base_stress.append(-1 * atoms.get_stress(voigt=False) * TO_KBAR)

        dft_energy.append(atoms.info['DFT_energy'] / len(atoms))
        dft_forces.append(atoms.arrays['DFT_forces'])
        dft_stress.append(-1 * atoms.info['DFT_stress'] * TO_KBAR)

    return {
        "dft": {
            "energy": dft_energy,
            "forces": flatten_arrays(dft_forces),
            "stress": flatten_arrays(dft_stress)
        },
        "base": {
            "energy": base_energy,
            "forces": flatten_arrays(base_forces),
            "stress": flatten_arrays(base_stress)
        },
        "fine": {
            "energy": ft_energy,
            "forces": flatten_arrays(ft_forces),
            "stress": flatten_arrays(ft_stress)
        }
    }

# ----------------------------- EOS Analysis -----------------------------
def atom_oneshot(atoms, calc):
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    stress = atoms.get_stress()
    calc_results = {"energy": energy, "forces": forces, "stress": stress}
    calculator = SinglePointCalculator(atoms, **calc_results)
    return calculator.get_atoms()

def atom_cell_relax(atoms, calc, logfile="-"):
    atoms.calc = calc
    cf = UnitCellFilter(atoms, hydrostatic_strain=True)
    opt = LBFGS(cf, logfile=logfile)
    opt.run(fmax=0.02, steps=1000)
    return atoms
def make_eos_structures(relaxed):
    relaxed_cell = relaxed.get_cell()
    relaxed_lat = relaxed_cell.lengths()[0]
    eos_structures = []
    for strain in np.linspace(-0.05, 0.05, 11):
        strained_lat = relaxed_lat * (1+strain)
        relaxed.set_cell([strained_lat]*3, scale_atoms=True)
        eos_structures.append(deepcopy(relaxed))
    return eos_structures

def get_eos_and_volume(eos_list):
    en_list = [atoms.get_potential_energy() for atoms in eos_list]
    vol_list = [atoms.get_volume() for atoms in eos_list]
    rel_en_list = np.array(en_list) - min(en_list)
    return rel_en_list, vol_list

# ----------------------------- Main Execution -----------------------------
def main():
    assert os.path.exists(DATA_PATH) and os.path.exists(os.path.join(DATA_PATH, 'evaluation'))

    fine_tuned_calc = SevenNetCalculator(os.path.join(WORKING_DIR, 'checkpoint_fine_tuned.pth'))
    base_calc = SevenNetCalculator('7net-0')

    traj = load_trajectory(os.path.join(DATA_PATH, 'evaluation/test_md.extxyz'))
    print(f"Loaded {len(traj)} frames from trajectory.")

    results = evaluate_parity(traj, fine_tuned_calc, base_calc)

    density_colored_scatter_plot(results["dft"]["energy"], results["base"]["energy"],
                                 results["dft"]["forces"], results["base"]["forces"],
                                 results["dft"]["stress"], results["base"]["stress"],
                                 '7net-0')

    density_colored_scatter_plot(results["dft"]["energy"], results["fine"]["energy"],
                                 results["dft"]["forces"], results["fine"]["forces"],
                                 results["dft"]["stress"], results["fine"]["stress"],
                                 'fine-tuned')

    os.makedirs('eos', exist_ok=True)
    atoms_list = read(os.path.join(DATA_PATH, 'evaluation/eos.extxyz'), ':')
    most_stable_idx = np.argmin([at.get_potential_energy() for at in atoms_list])
    atoms = atoms_list[most_stable_idx]

    print("Relaxing with fine-tuned potential...")
    ft_relaxed = atom_cell_relax(deepcopy(atoms), fine_tuned_calc, './eos/ft_relax_log.txt')

    print("Relaxing with base potential...")
    base_relaxed = atom_cell_relax(deepcopy(atoms), base_calc, './eos/base_relax_log.txt')

    ft_eos_structures = make_eos_structures(ft_relaxed)
    base_eos_structures = make_eos_structures(base_relaxed)

    ft_eos_oneshot = [atom_oneshot(stct, fine_tuned_calc) for stct in ft_eos_structures]
    base_eos_oneshot = [atom_oneshot(stct, base_calc) for stct in base_eos_structures]

    write('./eos/ft_eos.extxyz', ft_eos_oneshot)
    write('./eos/base_eos.extxyz', base_eos_oneshot)

    dft_eos, dft_vol = get_eos_and_volume(read(os.path.join(DATA_PATH, 'evaluation/eos.extxyz'), ':'))
    ft_eos, ft_vol = get_eos_and_volume(read('./eos/ft_eos.extxyz', ':'))
    base_eos, base_vol = get_eos_and_volume(read('./eos/base_eos.extxyz', ':'))

    plt.figure(figsize=(10/2.54, 8/2.54))
    plt.plot(dft_vol, dft_eos, label='DFT')
    plt.plot(ft_vol, ft_eos, label='Fine-tuned')
    plt.plot(base_vol, base_eos, label='SevenNet-0')
    plt.xlabel(r"Volume ($\\rm{\\AA}^3$)")
    plt.ylabel("Relative energy (eV)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("parity_plot_DFT_ML-eos.png", bbox_inches="tight", dpi=300)
    plt.show()

if __name__ == "__main__":
    main()

