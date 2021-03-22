# G. Stacking fault energy, Screw dislocations 1/2<111> for bcc High Entropy Alloys [VASP, ATAT MCSQS code]

[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
![versions](https://img.shields.io/pypi/pyversions/Django?color=green&label=python&style=plastic)
[![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html)

**[1] Generate Stacking fault energy & Screw dislocation for bcc crystals**

--> b is a burger vector (b = a/2[111])
GSFE = (E_fault - E_perfect)/Area

**[2] Local-lattice-distortion HEA Alloys**

Analysis of atomic mismatch for High Entropy Alloys. Various definition exists but I have chosen the one given in the paper referenced in the script. These definitions are arbitrary.

The script reads VASP POSCAR & CONTCAR file for initial and final coordinates and then analyse the ions drift from its initial position and compute the atomic mismatch.

**[3] Generating Random Structure using SQS technique or Short range order (SRO) method**

Python script to generate BCC/FCC/HCP random structures using SQS technique. This program is the MODIFICATION of the NANOHUB code.

https://icet.materialsmodeling.org/advanced_topics/sqs_generation.html 

https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/

**[4] Code to convert file generated from MCSQS ATAT code to VASP POSCAR file.**

**[5] Monte Carlo code for generating strucutre with lower Binding energy.**

**[6] https://www.ctcms.nist.gov/potentials/atomman/tutorial/04.3._Dislocation_analysis_tools.html**

**CITATION toward this work should be acknowledged in the publication.**
**USE AT YOUR OWN RISK. NOT EVEN IMPLIED WARRANTY WHATSOEVER**

NB:: For generating the structure you will need to install atomsk. The link is:
https://github.com/pierrehirel/atomsk/ 
it generates POSCAR files and they can be run independently with VASP and after relaxation in the z direction
the energy can be calulated and GSFE can be plotted against normalized burger vector.

**For Elastic constant**
http://wolf.ifj.edu.pl/elastic/index.html


https://shields.io/
