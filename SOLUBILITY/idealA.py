#!/usr/bin/env python

import numpy as np
import random, sys
from fortran_lib import solid_monte_carlo

#### Fundamental Variables ####
Temperature = 150 # kelvin
kspring = 6666.67 # kT/angstrom
gro_file = "solid.gro"
#rho = 0.01878
#iseed = random.randint(0,1000000)
# Get number of atoms
nAtom = 0
try:
    f = open(gro_file,'rb')
    lines = f.readlines()
    f.close()
except:
    raise IOError('No coordinate file')
nTot = lines[1] # line in the gro file that contains the total number of atoms
# Get all positions in the crystal lattice
AllPos = np.genfromtxt(gro_file, dtype=None, skip_header=2, skip_footer=1)
# Save coordinates of a single molecule for Monte Carlo calculations
Coord = []
for elem in AllPos:
    if '1MOL' in elem: # Substitute for '1FIX' in other inputs
        Coord.append(np.array([elem[3], elem[4], elem[5]]))
Coord = np.array(Coord)
nAtom = len(Coord)
nMol = float(nTot) / nAtom # number of molecules
centralize = np.tile(Coord[0],(len(Coord),1))
Pos = Coord - centralize # Rotations will be done about the origin. This ensures the fixed atom is there.
Pos = 10 * Pos # Converts to angstrom

print(Pos)
print(nAtom,nMol)

A0, error, I1, I2, dmax, cmin = solid_monte_carlo(Pos, kspring, nMol, rho=0.01878)

output = """

Helmholtz free energy of the ideal Einstein Molecule
A0 = %.5f kT at %.2f K
log(I1) = %.8f
log(I2) = %.8f

Maximum distance a molecule can move in the lattice:
dmax = %.5f angstrom

Maximum value of theta that contributes:
theta = %.5f rad 


""" % (A0, Temperature, I1, I2, dmax, np.cos(cmin))

out = open("resultA0.dat","w+")
out.write( output )
out.close()
print( output )
