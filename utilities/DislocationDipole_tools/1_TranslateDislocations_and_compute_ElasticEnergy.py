#!/usr/bin/env python3

"""
############################################################################
# Author: Asif Iqbal
# Date created: 2020-12-13
# Last modified: 2024-05-07
# Description: This script calculates the elastic energy of screw dislocations
#              by sampling different local chemical environments within a supercell.
#              The atoms are translated while keeping the screw dislocation fixed.
# Note: The POSCAR file should be in Cartesian coordinates.
############################################################################
"""

import numpy as np
import os, sys, re
import shutil
from tqdm import tqdm
from glob import glob

os.system("rm -rf dis_* iniTMP_* INPUT_GEO_CORE")

def read_POSCAR():
    pos = []
    kk = []
    lattice = []
    sum = 0
    with open("POSCAR_perfect_wo_dislo", "r") as file:
        firstline = file.readline()  # IGNORE first line comment
        alat = float(file.readline())  # scale
        
        Latvec1 = file.readline().split()
        Latvec2 = file.readline().split()
        Latvec3 = file.readline().split()
        
        elementtype = file.readline()
        atomtypes = file.readline()
        Coordtype = file.readline().split()
        
        if Coordtype[0] == "Direct" or Coordtype[0] == "direct":
            os.exit("First, Convert to Cartesian!")
            
        #nat = [int(i) for i in atomtypes.split()]
        nat = list(map(int, atomtypes.split()))
        n_atoms = np.sum(nat)    
            
        for x in range(int(n_atoms)):
            coord = [float(i) for i in file.readline().split()]
            pos = pos + [coord]
    
    Latvec1 = list(map(float, Latvec1))
    Latvec2 = list(map(float, Latvec2))
    Latvec3 = list(map(float, Latvec3))
    print(Latvec1)
    print(Latvec2)
    print(Latvec3)
    
    return (
        n_atoms,
        pos,
        firstline,
        alat,
        Latvec1,
        Latvec2,
        Latvec3,
        elementtype,
        atomtypes,
        Coordtype[0],
    )


def write_result(
    i,
    j,
    cnt,
    firstline,
    alat,
    Latvec1,
    Latvec2,
    Latvec3,
    elementtype,
    atomtypes,
    Coordtype,
    pos,
    n_atoms):
    
    # Create the file and write data
    filename = f"iniTMP_{cnt:02d}"
    
    with open(filename, "w") as fdat1:
        fdat1.write(f"{firstline}")  # Comment line in POSCAR
        fdat1.write(f"{alat:.5f}\n")
        fdat1.write(f"{Latvec1[0]:.11f} {Latvec1[1]:.11f} {Latvec1[2]:.11f}\n")
        fdat1.write(f"{Latvec2[0]:.11f} {Latvec2[1]:.11f} {Latvec2[2]:.11f}\n")
        fdat1.write(f"{Latvec3[0]:.11f} {Latvec3[1]:.11f} {Latvec3[2]:.11f}\n")
        fdat1.write(f"{elementtype}")
        fdat1.write(f"{atomtypes}")
        fdat1.write(f"{Coordtype}\n")

        # Displace the cell in the "X" & "Y" direction.
        for x in range(n_atoms):
            fdat1.write(f"{pos[x][0] + i:.12f} {pos[x][1] + j:.12f} {pos[x][2]:.12f}\n")


# MAIN
(
    n_atoms,
    pos,
    firstline,
    alat,
    Latvec1,
    Latvec2,
    Latvec3,
    elementtype,
    atomtypes,
    Coordtype,
) = read_POSCAR()

Ax = float(Latvec1[0])
Ay = float(Latvec2[1])
Az = float(Latvec3[2])

cnt = 0
cx, cy = 0, 0

b = (3.404 / 2) * np.linalg.norm([1, 1, 1]) # b = a0/2 |<111>|
a0 = (2 * b) / np.sqrt(3) # 2*b / |<111>|

print(f"SYSTEM={elementtype.split()}, \n#atoms={n_atoms}, \nb={b}")
print(f"a0 = {a0:6.4f}")

os.makedirs("INPUT_GEO_CORE", exist_ok=True)

# === Tunable parameters 
with open("FILE_INFO.csv", "w") as ff:
    ff.write(f"Ax[A],Ay[A],Elastic_energy[eV/length],counter\n")
    batch_size = 2
    i_values = np.linspace(0, Ax, num=batch_size)
    j_values = np.linspace(0, Ay, num=batch_size)
    print(f"# of point => {len(i_values)*len(j_values)}")
    
    pattern = r'Total elastic energy :.*=.*='
    
    for i in tqdm(i_values):
        for j in j_values:
            write_result(i,j,cnt,firstline,alat,Latvec1,Latvec2,Latvec3,elementtype,atomtypes,Coordtype,pos,n_atoms)

            L = f"dis_X{cx:02d}_Y{cy:02d}_{cnt:02d}"
            K = f"disFIN_X{cx:02d}_Y{cy:02d}_{cnt:02d}"
            
            shutil.rmtree(L, ignore_errors=True)  # overwrite a directory
            os.makedirs(L, exist_ok=True)

            shutil.copy(f"iniTMP_{cnt:02d}", L)
            shutil.copy("INCAR", os.path.join(L, "INCAR"))
            shutil.copy("POTCAR", os.path.join(L, "POTCAR"))
            shutil.copy("KPOINTS", os.path.join(L, "KPOINTS"))
            shutil.copy("input_dislo.babel", os.path.join(L, "input_dislo.babel"))
            
            os.chdir(L)

            shutil.copy(f"iniTMP_{cnt:02d}", "POSCAR_perfect")
            
            # return from babel POSCAR
            os.system("babel input_dislo.babel > ELASTIC_ENE.dat")
            
            with open('ELASTIC_ENE.dat', 'r') as CE:
                for line in CE:
                    if re.search(pattern, line):
                        match = re.search(r'= *(\S+) *=', line)
                        if match:
                            elastic_energy = float(match.group(1))/Az # eV per unit length
                            break  

            shutil.copy("POSCAR", K)
            os.chdir("../")
            
            shutil.copy(os.path.join(L, K), "INPUT_GEO_CORE")
            
            ff.write(f"{i:12.6f},{j:12.6f},{elastic_energy:12.6f},{cnt:3d}\n")

            cnt += 1
            cy += 1
        cy = 0  # reset counter
        cx += 1
        
for dd in glob('iniTMP_*'):        
    shutil.move(dd, "INPUT_GEO_CORE")

print("DISPLACED FILES HAS BEEN GENERATED in the X=<112>, Y=<110>, Z=<111>")


