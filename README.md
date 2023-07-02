Scripts to analyze Screw dislocations, GSFE, Frenkel Kontorova, NEB, and SRO for High Entropy Alloys (HEA)
ML and DNN and linear regression scripts
ML-based potential for HEAs (ternary to quinary) is in progress

[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
![versions](https://img.shields.io/pypi/pyversions/Django?color=green&label=python&style=plastic)
[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html)

**USE AT YOUR OWN RISK. NOT EVEN AN IMPLIED WARRANTY, WHATSOEVER!

Scripts that are included in this package:
_________________
1. Generate Stacking fault energy & Screw dislocation for bcc crystals

   --> b is a burger vector (b = a/2[111]) and GSFE is calculated as
   GSFE = (E_fault - E_perfect)/Area

   NB:: For generating the structure you will need to install atomsk. The link is:
https://github.com/pierrehirel/atomsk/ 
it generates POSCAR files and they can be run independently with VASP and after relaxation in the z-direction,
the energy can be calculated. GSFE can be plotted against a normalized burger vector.

2. Local-lattice-distortion HEA Alloys

   Analysis of atomic mismatch for High Entropy Alloys. Various definitions exist, but I have chosen the one in the paper referenced in the script. These definitions are arbitrary. The script reads VASP POSCAR & CONTCAR file for initial and final coordinates and then analyzes the ions drift from its initial position and compute the atomic mismatch.

3. Generating Random Structure using SQS technique method: Python script to generate BCC/FCC/HCP random structures. This program is the MODIFICATION of the NANOHUB code.

4. Code to convert file generated from MCSQS ATAT code to VASP POSCAR file

5. Monte Carlo code for generating structure with lower Binding energy

6. For Dislocation analysis: 

https://www.ctcms.nist.gov/potentials/atomman/tutorial/04.3._Dislocation_analysis_tools.html

7. Creation of Screw dislocations
8. Finding SRO with pyscal code
9. Find the core energy of screw dislocation 

**BIBLIOGRAPHY
___________________________
1. https://icet.materialsmodeling.org/advanced_topics/sqs_generation.html; 
2. https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/; 
3. https://aip.scitation.org/doi/10.1063/5.0014094; 
4. https://pyscal.org/en/latest/
5. http://emmanuel.clouet.free.fr/presentation.html
6. https://github.com/atomistic-machine-learning/schnetpack
7. https://github.com/libAtoms/QUIP
8. https://github.com/ACEsuit/mace
9. https://github.com/Liu-group/AutoSolvate
10. https://github.com/jbuckeridge/cplap
11. https://github.com/MaterSim/PyXtal_FF
12. https://github.com/mir-group/nequip
13. https://github.com/uw-cmg/MAST-ML



https://shields.io/

#VASP, #ATAT, #MCSQS, #ASE, #Pymatgen, #Python, #DNN, #Keras, #TensorFlow
