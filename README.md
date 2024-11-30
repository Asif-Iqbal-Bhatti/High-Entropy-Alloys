# <span style="color: green">Selected Python Codes for High Entropy Alloys Analysis</span>

[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)  
![versions](https://img.shields.io/pypi/pyversions/Django?color=green&label=python&style=plastic)  
[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html)  

### DISCLAIMER: Please be sure to use it at your risk. THERE IS NO IMPLIED WARRANTY WHATSOEVER!

## Scripts used in this paper are uploaded to this repository https://doi.org/10.1016/j.commatsci.2024.113196

## Overview
_> This GitHub repository contains a collection of Python scripts for analyzing various properties of High Entropy Alloys (HEAs). The scripts cover a range of topics, including screw dislocations, Generalized Stacking Fault Energy (GSFE), Frenkel Kontorova, Nudged Elastic Band (NEB), and Short-Range Order (SRO). Additionally, it includes Machine Learning (ML), Deep Neural Networks (DNN), and linear regression codes to study HEAs. Furthermore, a work-in-progress ML-based potential for HEAs, from ternary to quinary compositions, is also included._

### Scripts Included
> Generate Stacking Fault Energy & Screw Dislocation for BCC Crystals:
This script calculates the GSFE and screw dislocation for bcc crystals using the formula:
> '''
GSFE = (E_fault - E_perfect) / Area
> '''
Atomsk needs to be installed to generate the crystal structures. The link to Atomsk is provided (https://github.com/pierrehirel/atomsk/). The generated POSCAR files can be independently run with VASP after relaxation in the z-direction to calculate the energy. GSFE can be plotted against a normalized burger vector.

### Local Lattice Distortion in HEA Alloys:
> This script analyzes atomic mismatches for high-entropy alloys. It reads the VASP POSCAR and CONTCAR files for initial and final coordinates, computes the atomic drift from its initial position, and calculates the atomic mismatch. The definitions used in the script are based on a referenced paper, although various definitions exist.
This script analyzes atomic mismatches for high-entropy alloys. It reads the VASP POSCAR and CONTCAR files for initial and final coordinates, computes the atomic drift from its initial position, and calculates the atomic mismatch. The definitions used in the script are based on a referenced paper, although various definitions exist.

### Generating Random Structure Using the SQS Technique:
> Python script to generate BCC/FCC/HCP random structures using the modified NANOHUB code.

Code to Convert File Generated from MCSQS ATAT Code to VASP POSCAR File:
This script converts files generated from the MCSQS ATAT code to VASP POSCAR format.

### Monte Carlo Code:
> Monte Carlo code for generating structures with lower binding energy

### Dislocation Analysis:
> Script for analyzing dislocations using tools available at https://www.ctcms.nist.gov/potentials/atomman/tutorial/04.3._Dislocation_analysis_tools.html.

### Creation of Screw Dislocations:
> Script for creating screw dislocations.

### Finding Short-Range Order (SRO) with Pyscal Code:
> Script using the pyscal code to find short-range Orders in the system

### Finding the Core Energy of Screw Dislocation:
> Script to determine the core energy of a screw dislocation.

## Bibliography
________
The scripts are based on concepts and techniques from the following references:

https://icet.materialsmodeling.org/advanced_topics/sqs_generation.html  
https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/  
https://aip.scitation.org/doi/10.1063/5.0014094  
https://pyscal.org/en/latest/  
http://emmanuel.clouet.free.fr/presentation.html  
https://github.com/atomistic-machine-learning/schnetpack  
https://github.com/libAtoms/QUIP  
https://github.com/ACEsuit/mace  
https://github.com/Liu-group/AutoSolvate  
https://github.com/jbuckeridge/cplap  
https://github.com/MaterSim/PyXtal_FF  
https://github.com/mir-group/nequip  
https://github.com/uw-cmg/MAST-ML  
http://ann.atomistic.net/links/  
https://openkim.org/  
https://github.com/materialsvirtuallab/matgl  
https://www.doitpoms.ac.uk/tlplib/dislocations/index.php
https://github.com/gcmt-group/sod

## Notes

_This repository is provided "as is," and no guarantees or warranties are associated with the code.
Users are advised to exercise caution and thoroughly review the scripts before using them.
The repository is actively maintained, and updates may be available to enhance or improve the scripts._

## Tags

#VASP #ATAT #MCSQS #ASE #Pymatgen #Python #DNN #Keras #TensorFlow

Feel free to use and contribute to this repository. For more information, kindly visit the project's GitHub page.
