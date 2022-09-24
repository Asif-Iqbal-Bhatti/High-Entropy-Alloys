#!/usr/bin/env python3

#===========================================
# CONVERT wt% to at%
# EXP:: Ti-25Nb-8Zr (in wt%)
# could be generalized to any composition
#===========================================

from ase.data import atomic_numbers, atomic_names, atomic_masses, covalent_radii

Avogadro = 6.022 * 10**23

Nb = atomic_masses[atomic_numbers['Nb']] # 93
Zr = atomic_masses[atomic_numbers['Zr']] # 91.2
Ti = atomic_masses[atomic_numbers['Ti']] # 48

Atomic_weight = [Nb, Zr, Ti]

nb = 25 * Avogadro / Atomic_weight[0]
zr = 8 * Avogadro / Atomic_weight[1]
ti = 67 * Avogadro / Atomic_weight[2]

nb_at = nb / (nb + zr + ti)
zr_at = zr / (nb + zr + ti)
ti_at = ti / (nb + zr + ti)

gg = {'Nb': nb_at, 'Zr': zr_at, 'Ti': ti_at}

for i, j in gg.items():
	print(f'Therefore, at% of {i} {j *100:3.2f}%')


