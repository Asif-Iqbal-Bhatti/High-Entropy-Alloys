#!/usr/bin/env python3

'''
#===========================================
# AUTHOR:: AIB_EM
# CONVERT wt% to at%
# EXP:: Ti-25Nb-8Zr (in wt%)
#===========================================
'''

import re, os
from ase.data import atomic_numbers, atomic_names, atomic_masses, covalent_radii
from ase.formula import Formula

Avogadro = 6.022 * 10**23

tot = 0
a, g = {}, {}
at_weight = {}

chem_for = 'Ti-25Nb-8Zr-5F'

e = re.findall('[A-Z][a-z]?', str(chem_for))
w = chem_for.split('-')

for i, j in enumerate(w):
	if list(Formula(j).count().values())[-1] > 1:
		a.update(Formula(j).count())
rest = list(set(e) - set(a))[0]

for k, v in a.items():
	tot += v
a[rest] = 100 - tot # Weight is out of 100

for k, _ in a.items():
	at_weight[k] = atomic_masses[atomic_numbers[k]]

print(f'WT% IN FORMULA -> {a}')
print(f'ATOMIC WT% -> {at_weight}')

for k, v in a.items():
	g[k] = (a[k] * Avogadro) / float(at_weight[k])
all = sum(list(g.values()))

for k, v in a.items():
	d = g[k] / all
	print(f"Therefore, at% of {k} is {d*100:3.3f}%")


#=========== VERIFICATION
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
	#print(f'Therefore, at% of {i} {j *100:3.2f}%')
	pass
