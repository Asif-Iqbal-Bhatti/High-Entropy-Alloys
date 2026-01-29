#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR :: AIBEM
DATED  :: 29/01/2026
DOCUMENTATION :: Fast PM–SRO (and WC–SRO) using ASE neighbor list, supporting K shells.
                 REWRITTEN version of PMSRO in PYSCAL2

Usage example
-------------
from PMSROCalculatorPyScale2 import PMSROCalculator

calc = PMSROCalculator(
    element_group=['Mn', 'Ti', 'Nb'],  # or None for all species
    delta=True                          # PM–SRO for i==j (de Fontaine); False -> pure WC–SRO
)

# Provide strictly increasing shell cutoffs [r1, r2, ..., rK]
calc.compute(shell_cutoffs=[3.2, 4.4, 5.6])    # 3 shells
print(calc.table.round(4))                     # tidy DataFrame (1st, 2nd, 3rd, Total)
print(calc.alpha_total)                        # Z-weighted total alpha_ij across K shells
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import Iterable, List, Optional, Tuple

from ase.io import read
from ase.neighborlist import neighbor_list
# Optional: faster alternative neighbor finder with identical signature:
# from matscipy.neighbours import neighbour_list as ms_neighbour_list    # [2](https://journals.iucr.org/a/issues/2012/05/00/ib5006/ib5006.pdf)

@dataclass
class PMSROCalculator:
    element_group: Optional[List[str]] = None
    delta: bool = True  # True -> PM–SRO for i==j (de Fontaine), False -> pure WC–SRO
    atoms: Optional[object] = field(default=None, init=False)
    elements: Optional[List[str]] = field(default=None, init=False)
    el_to_idx: Optional[dict] = field(default=None, init=False)
    conc: Optional[np.ndarray] = field(default=None, init=False)      # c_j within group
    shell_cutoffs: Optional[List[float]] = field(default=None, init=False)
    alpha_shells: Optional[List[np.ndarray]] = field(default=None, init=False)  # [alpha^(1), ..., alpha^(K)], each (G,G)
    alpha_total: Optional[np.ndarray] = field(default=None, init=False)         # (G,G)
    table: Optional[pd.DataFrame] = field(default=None, init=False)


    def from_poscar(self, path: str) -> None:
        """Load structure from a POSCAR/CONTCAR path (read by ASE)."""
        self.atoms = read(path)

    def from_atoms(self, atoms) -> None:
        """Load structure from an ASE Atoms object."""
        self.atoms = atoms

    # ----------------------------
    # Main compute
    # ----------------------------
    def compute(self, shell_cutoffs: Iterable[float]) -> None:
        """
        Compute PM–SRO/WC–SRO over K shells defined by monotonically increasing cutoffs.

        Parameters
        ----------
        shell_cutoffs
            Iterable of strictly increasing radii [r1, r2, ..., rK]. Shell s is (r_{s-1}, r_s] with r_0 := 0.
        """
        if self.atoms is None:
            raise ValueError("No structure loaded. Use from_poscar() or from_atoms().")

        # Normalize/validate cutoffs
        cuts = np.array(list(shell_cutoffs), dtype=float).ravel()
        if cuts.ndim != 1 or len(cuts) < 1 or not np.all(np.diff(cuts) > 0):
            raise ValueError("shell_cutoffs must be a 1D strictly increasing sequence of radii (length >= 1).")
        self.shell_cutoffs = cuts.tolist()
        K = len(cuts)

        # Symbols and group selection
        symbols_all = np.array(self.atoms.get_chemical_symbols())
        if self.element_group is None:
            #self.elements = sorted(set(symbols_all.tolist()))
            self.elements = list(dict.fromkeys(symbols_all))
        else:
            # keep only elements that actually exist in structure
            present = sorted(set(s for s in self.element_group if (symbols_all == s).any()))
            if not present:
                raise ValueError("Selected element_group not present in structure.")
            self.elements = present
        G = len(self.elements)
        self.el_to_idx = {el: i for i, el in enumerate(self.elements)}

        # Mask atoms to the selected group (e.g., TM-only for DRX)
        group_mask = np.isin(symbols_all, self.elements)
        idx_map = np.where(group_mask)[0]
        if len(idx_map) == 0:
            raise ValueError("No atoms found for the selected element_group.")

        # Concentrations within the selected group
        group_symbols = symbols_all[idx_map]
        group_counts = np.array([np.count_nonzero(group_symbols == el) for el in self.elements], dtype=float)
        self.conc = group_counts / group_counts.sum()     # c_j

        # Build neighbor list once up to the outermost cutoff rK (directed pairs i->j)
        r_max = cuts[-1]
        # i_all, j_all, d_all = ms_neighbour_list('ijd', self.atoms, r_max)  # matscipy alternate (often very fast)  [2](https://journals.iucr.org/a/issues/2012/05/00/ib5006/ib5006.pdf)
        i_all, j_all, d_all = neighbor_list('ijd', self.atoms, r_max)        # ASE neighbor list  [1](https://pubs.rsc.org/en/content/articlelanding/2023/ta/d3ta02475j)

        # Keep pairs where both center and neighbor are in the chosen group
        mask_group = np.isin(i_all, idx_map) & np.isin(j_all, idx_map)
        i_all, j_all, d_all = i_all[mask_group], j_all[mask_group], d_all[mask_group]

        # Remap absolute indices to compact group positions (0..Ng-1) then to species indices (0..G-1)
        pos_in_group = -np.ones(len(symbols_all), dtype=int)
        pos_in_group[idx_map] = np.arange(len(idx_map))
        gi = pos_in_group[i_all]
        gj = pos_in_group[j_all]
        ti = np.fromiter((self.el_to_idx[symbols_all[idx_map[g]]] for g in gi), count=len(gi), dtype=int)  # type of i
        tj = np.fromiter((self.el_to_idx[symbols_all[idx_map[g]]] for g in gj), count=len(gj), dtype=int)  # type of j

        # Build K shell masks: (0, r1], (r1, r2], ..., (r_{K-1}, rK]
        prev = 0.0
        shell_masks: List[np.ndarray] = []
        for r in cuts:
            m = (d_all > prev) & (d_all <= r)
            shell_masks.append(m)
            prev = r

        # Per-shell counts N_ij and row sums Z_i
        Nijs: List[np.ndarray] = []
        Zis:  List[np.ndarray] = []
        for m in shell_masks:
            if not np.any(m):
                N_ij = np.zeros((G, G), dtype=float)
                Zi   = np.zeros(G, dtype=float)
            else:
                # Count directed pairs by species (vectorized)
                N_ij = np.zeros((G, G), dtype=float)
                np.add.at(N_ij, (ti[m], tj[m]), 1.0)
                Zi = N_ij.sum(axis=1)
            Nijs.append(N_ij)
            Zis.append(Zi)

        # Compute alpha_ij for each shell
        alpha_shells: List[np.ndarray] = []
        c = self.conc
        for N_ij, Zi in zip(Nijs, Zis):
            with np.errstate(divide='ignore', invalid='ignore'):
                p_ij = N_ij / Zi[:, None]  # rows with Zi==0 -> NaN
            # start from WC–SRO for all pairs
            alpha = 1.0 - (p_ij / c[None, :])
            if self.delta:
                # PM–SRO correction for i==j: alpha_ii = (p_ii - c_i)/(1 - c_i)
                for k in range(G):
                    if np.isfinite(p_ij[k, k]):
                        alpha[k, k] = (p_ij[k, k] - c[k]) / (1.0 - c[k])
            alpha_shells.append(alpha)

        self.alpha_shells = alpha_shells

        # Z-weighted total across shells, row-wise
        alpha_total = np.zeros((G, G), dtype=float)
        for irow in range(G):
            weights = np.array([Zi[irow] for Zi in Zis], dtype=float)
            good = (weights > 0) & np.isfinite(weights)
            if not np.any(good):
                alpha_total[irow, :] = np.nan
            else:
                num = sum(w * np.nan_to_num(a[irow, :], nan=0.0) for w, a in zip(weights, alpha_shells))
                den = weights[good].sum()
                alpha_total[irow, :] = num / den
        self.alpha_total = alpha_total

        # Build tidy table with K shells + Total
        rows = []
        idx = []
        for i, ei in enumerate(self.elements):
            for j, ej in enumerate(self.elements):
                idx.append(f"{ei}-{ej}")
                vals = [a[i, j] for a in alpha_shells]   # shell 1..K
                # Z-weighted total equals self.alpha_total[i, j]
                rows.append(vals + [alpha_total[i, j]])
        colnames = [f"{s+1}st" if s == 0 else (f"{s+1}nd" if s == 1 else (f"{s+1}rd" if s == 2 else f"{s+1}th"))
                    for s in range(K)]
        colnames.append("Total")
        self.table = pd.DataFrame(rows, index=idx, columns=colnames)

    # ----------------------------
    # Convenience getters
    # ----------------------------
    def get_shell_matrix(self, shell_index: int) -> np.ndarray:
        """Return alpha^(shell_index) matrix (0-based index)."""
        if self.alpha_shells is None:
            raise ValueError("Nothing computed yet. Call compute().")
        return self.alpha_shells[shell_index]

    def get_total_matrix(self) -> np.ndarray:
        if self.alpha_total is None:
            raise ValueError("Nothing computed yet. Call compute().")
        return self.alpha_total