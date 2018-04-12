#!/usr/bin/env python

# the file passed to the program is the gaussian output file from which it reads
# the geometry of the molecule
# assumes that inertia.dat is present in the current directory and reads
# the inertia data

import sys
import readGeomFc
import os
from numpy import *
import geomUtility

# initialize parameters
#seed = 200
numIter = 20000
theory = ' b3lyp/6-31G(d) '
numprocessors = 8
memory = ' 1000MB '
multi = ' 2'

inertia = open('inertia.dat', 'r')
result = open('dihed_energy.out', 'w')

outputName = '0.log'
energy = readGeomFc.readHFEnergy(outputName)
outFile = open('0.log', 'r')
(geom, Mass) = readGeomFc.readGeom(outFile)

(rotors) = readGeomFc.readGeneralInertia(inertia, Mass)
numRotors = len(rotors) - 1


atoms = readGeomFc.getAtoms(Mass)
diheds = numRotors * [0]
for i in range(numRotors):
    result.write('%10.3f' % diheds[i])
result.write('%14.7f' % energy + '\n')


for i in range(12**(numRotors)):
    irem = i
    for j in range(len(rotors) - 1):
        diheds[j] = (irem - irem / 12 * 12) * 30.0
        irem = irem / 12
        result.write('%10.3f' % diheds[j])
    result.write('\n')

    newGeom = geomUtility.rotateDihedrals(geom, rotors, diheds, Mass)

    file = open(str(i + 1) + '.com', 'w')
    file.write('%mem=' + memory + '\n')
    file.write('%nproc=' + str(numprocessors) + '\n')
    file.write('#' + theory + ' nosym\n\nTitle Card Required\n\n0' + multi + '\n')
    for j in (range(len(Mass))):
        file.write(atoms[j] + '  ' + str(float(newGeom[j, 0])) + '  ' +
                   str(float(newGeom[j, 1])) + '  ' + str(float(newGeom[j, 2])) + '\n')
    file.write('\n\n')
    file.close()
result.close()
