#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from ase.data import atomic_numbers
import pandas as pd

# ----- inputs -----
base = sys.argv[1] if len(sys.argv) > 1 else "."
elements = ["Cl", "Li", "P", "S"]

# MACE-MP reference energy
ref_energy_dict = {
    1: -1.11734008, 3: -0.29754725, 4: -0.01781697, 5: -0.26885011, 6: -1.26173507,
    7: -3.12438806, 8: -1.54838784, 9: -0.51882044, 11: -0.22883163, 12: -0.00951015,
    13: -0.21630193, 14: -0.8263903, 15: -1.88816619, 16: -0.89160769, 17: -0.25828273,
    19: -0.22697913, 20: -0.0927795, 21: -2.11396364, 22: -2.50054871, 23: -3.70477179,
    24: -5.60261985, 25: -5.32541181, 26: -3.52004933, 27: -1.93555024, 28: -0.9351969,
    29: -0.60025846, 30: -0.1651332, 31: -0.32990651, 32: -0.77971828, 33: -1.68367812,
    34: -0.76941032, 35: -0.22213843, 36: -0.0335879, 37: -0.1881724, 38: -0.06826294,
    39: -2.17084228, 40: -2.28579303, 41: -3.13429753, 42: -4.60211419, 43: -3.45201492,
    44: -2.38073513, 45: -1.46855515, 46: -1.4773126, 47: -0.33954585, 48: -0.16843877,
    49: -0.35470981, 50: -0.83642657, 51: -1.41101987, 52: -0.65740879, 53: -0.18964571,
    54: -0.00857582, 55: -0.13771876, 56: -0.03457659, 57: -0.45580806, 58: -1.3309175,
    59: -0.29671824, 60: -0.30391193, 61: -0.30898427, 62: -0.25470891, 63: -8.38001538,
    64: -10.38896525, 65: -0.3059505, 66: -0.30676216, 67: -0.30874667, 68: -0.31610927,
    69: -0.25190039, 70: -0.06431414, 71: -0.31997586, 72: -3.52770927, 73: -3.54492209,
    74: -4.65658356, 75: -4.70108713, 76: -2.88257209, 77: -1.46779304, 78: -0.50269936,
    79: -0.28801193, 80: -0.12454674, 81: -0.31737194, 82: -0.77644932, 83: -1.32627283,
    89: -0.26827152, 90: -0.90817426, 91: -2.47653193, 92: -4.90438537, 93: -7.63378961,
    94: -10.77237713
}

# ----- collect energies -----
vacuum_dirs = [d for d in os.listdir(base) if re.search(r"vac", d, re.IGNORECASE)]
vacuum_dirs.sort(key=lambda x: float(re.search(r"([0-9]+(?:\\.[0-9]+)?)", x).group(1)))

data = {el: [] for el in elements}
vac_values = []

for vdir in vacuum_dirs:
    vac = float(re.search(r"([0-9]+(?:\\.[0-9]+)?)", vdir).group(1))
    vac_values.append(vac)
    for el in elements:
        xml = os.path.join(base, vdir, el, "vasprun.xml")
        atoms = read(xml, index=-1)  # last ionic step
        e_calc = atoms.get_potential_energy()
        data[el].append(e_calc if e_calc is not None else float("nan"))


# ----- build a pandas table: Element, Z, Vacuum, Calc_E, Ref_E, Delta_E -----
rows = []
for el in elements:
    Z = atomic_numbers[el]
    e_ref = ref_energy_dict.get(Z, np.nan)
    for vac, e_calc in zip(vac_values, data[el]):
        dE = e_calc - e_ref if (np.isfinite(e_ref) and np.isfinite(e_calc)) else np.nan
        rows.append({
            "Element": el,
            "Z": Z,
            "Vacuum (Å)": vac,
            "Calc_E (eV)": e_calc,
            "Ref_E (eV)": e_ref,
            "ΔE = Calc - Ref (eV)": dE
        })

df = pd.DataFrame(rows)
# nice formatting for console
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.width", 120)
pd.set_option("display.precision", 6)

print(df.sort_values(["Element", "Vacuum (Å)"]).to_string(index=False))


# ----- plots -----
# Combined plot (log x-axis)
plt.figure(figsize=(8, 6), dpi=150)
for el in elements:
    plt.plot(vac_values, data[el], marker="o", label=el)
plt.xscale("log")
plt.xlabel("Vacuum (Å)")
plt.ylabel("Energy (eV)")
plt.title("Energy vs Vacuum")
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.legend()
plt.tight_layout()
plt.savefig("energies.png", bbox_inches="tight", dpi=300)

# Per-element subplots with tight y-limits
fig, axes = plt.subplots(2, 2, figsize=(10, 8), dpi=150)
axes = axes.flatten()
for i, el in enumerate(elements):
    ax = axes[i]
    y = np.array(data[el], dtype=float)
    ax.plot(vac_values, y, marker="o", linewidth=1.8)
    ax.set_xscale("log")
    ax.set_title(el)
    ax.set_xlabel("Vacuum (Å) [log]")
    ax.set_ylabel("Energy (eV)")
    ymin, ymax = np.nanmin(y), np.nanmax(y)
    margin = max(1e-3, 0.02 * max(1e-6, abs(ymax - ymin)))
    ax.set_ylim(ymin - margin, ymax + margin)
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
plt.suptitle("Energy vs Vacuum (Per-element convergence)", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("energies_convergence_subplots.png", bbox_inches="tight", dpi=300)

# Optional ΔE plot relative to reference
plt.figure(figsize=(8, 6), dpi=150)
for el in elements:
    Z = atomic_numbers[el]
    e_ref = ref_energy_dict.get(Z, np.nan)
    dy = np.array(data[el], dtype=float) - e_ref
    plt.plot(vac_values, dy, marker="o", linewidth=1.8, label=el)
plt.xscale("log")
plt.xlabel("Vacuum (Å) [log]")
plt.ylabel("ΔE = E_calc - E_ref (eV)")
plt.title("Convergence relative to reference energy")
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.legend()
plt.tight_layout()
plt.savefig("energies_delta_log.png", bbox_inches="tight", dpi=300)



