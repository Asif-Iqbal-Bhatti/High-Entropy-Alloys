# Stacking fault energy & Screw dislocations for bcc crystals for VASP setup // ATAT MCSQS code 

**[1] Generate Stacking fault energy & Screw dislocation for bcc crystals for VASP setup**

|||--> b means burger vector
|||--> b = a*1/2[111] ;|b| = axsqrt(3)/2
|||--> On the (110) plane the slip occur at <111> 

NB:: For executation you will need to install atomsk. The link is:
https://github.com/pierrehirel/atomsk/ 
it generates POSCAR files and these can be run independently with VASP and after relaxation in the z direction
the energy can be calulated and GSFE can be calculated and plotted against normalized burger vector.

GSFE = (E_fault - E_perfect)/Area

```
Some bcc materials (e.g. α-Fe) can contain up to 48 slip systems. 
There are six slip planes of type {110}, each with two <111> directions (12 systems). 
There are 24 {123} and 12 {112} planes each with one <111> direction (36 systems, 
for a total of 48). While the {123} and {112} planes are not exactly identical in 
activation energy to {110}, they are so close in energy that for all intents and 
purposes they can be treated as identical. example: 
slip plane and direction are (110) and [-111], respectively
```

**[2] Local-lattice-distortion HEA Alloys**

Analysis of atomic mismatch for High Entropy Alloys. Various definition exists but I have chosen the one given in the paper referenced in the script. These definitions are arbitrary.

The script reads VASP POSCAR & CONTCAR file for initial and final coordinates and then analyse the ions drift from its initial position and compute the atomic mismatch.

**[3] Generating Random Structure by using SQS technique**

Python script to generate BCC/FCC/HCP random structures using SQS technique. This program is the MODIFICATION of the NANOHUB code.
https://icet.materialsmodeling.org/advanced_topics/sqs_generation.html 

https://github.com/dgehringer/sqsgenerator 

https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/

Convert file generated from MCSQS ATAT code to VASP POSCAR file. 
**[4] Monte Carlo algorithm for generating strucutre with lower Binding energy

**CITATION toward this work should be acknowledged in the publication**
