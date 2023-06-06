#!/usr/bin/env python3

# USING ASE RELAXING POSITIONS AND LATTICE PARAMETERS
# AUTHOR:: Asif Iqbal - > AIB_EM
# https://wiki.fysik.dtu.dk/ase/ase/constraints.html#the-expcellfilter-class
# 

import os, sys
os.environ["OMP_NUM_THREADS"] = "8"
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
os.environ['TF_ENABLE_ONEDNN_OPTS'] = '1'
from ase.constraints import ExpCellFilter
from ase.optimize import LBFGS
from ase import Atoms
from ase.io import read
from m3gnet.models import M3GNet, M3GNetCalculator, Potential, Relaxer

myM3GNet = M3GNet.load('MP-2021.2.8-EFS')
myPotential = Potential(myM3GNet)
pot = M3GNetCalculator(
        potential=myPotential, 
        compute_stress = True, 
        stress_weight=0.01
)

atoms = read("POSCAR")
print(atoms)

atoms.calc = pot
ucf = ExpCellFilter(atoms)
opt = LBFGS(ucf)
opt.run(fmax=0.005)

atoms.write('CONTCAR')
