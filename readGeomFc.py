#!/usr/bin/env python
# takes the file as an argument and then reads the geometry and lower triangular part
# of the force constant matrix and the Atomic mass vector
'''
    Copyright (C) 2018, Sandeep Sharma and James E. T. Smith

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import os
from numpy import *
from Rotor import *
import pdb
import CanTherm
from Molecule import *
import re
import shutil
import getopt


def readInputFile(file, data):

    # read calculation type
    line = readMeaningfulLine(file)
    if line.split()[0].upper() == 'THERMO':
        data.CalcType = 'Thermo'
    elif line.split()[0].upper() == 'REAC':
        data.CalcType = 'Reac'
    else:
        print('first line of the input file is neither Reac or Thermo')
        exit()

# read reaction type
    line = readMeaningfulLine(file)
    if data.CalcType == 'Reac':
        if line.split()[0].upper() == 'BIMOL':
            data.ReacType = 'Bimol'
        elif line.split()[0].upper() == 'UNIMOL':
            data.ReacType = 'Unimol'
        else:
            print('ReactionType is not specified in the. Either Unimol or Bimol should be specified after keyword Read')
            exit()
        line = readMeaningfulLine(file)

# read Temperature range
    if line.split()[0].upper() == 'TLIST':
        line = readMeaningfulLine(file)
        numTemp = int(line.split()[0])
        i = 0
        while i < numTemp:
            line = readMeaningfulLine(file)
            tokens = line.split()
            i = i + len(tokens)
            for j in tokens:
                data.Temp.append(float(j))
        if len(data.Temp) > numTemp:
            print('More Temperatures than ', numTemp, ' are specified')

    elif line.split()[0].upper() == 'TRANGE':
        line = readMeaningfulLine(file)
        tokens = line.split()
        T0 = float(tokens[0])
        dT = float(tokens[1])
        numTemp = int(tokens[2])
        for i in range(numTemp):
            data.Temp.append(T0 + i * dT)
    else:
        print('Temperaure information not given either use keyword Tlist or Trange')
        exit()

# read scaling factor
    line = readMeaningfulLine(file)
    if (line.split()[0].upper() != 'SCALE'):
        print('Give a scaling factor')
    data.scale = float(line.split()[1])

# read molecular data
    if (data.CalcType == 'Thermo'):
        numMol = 1
    elif data.CalcType == 'Reac' and data.ReacType == 'Unimol':
        numMol = 2
    elif data.CalcType == 'Reac' and data.ReacType == 'Bimol':
        numMol = 3

    for i in range(numMol):
        if data.ReacType == 'Unimol' and i == 1 or data.ReacType == 'Bimol' and i == 2:
            molecule = Molecule(file, True)
        else:
            molecule = Molecule(file, False)
        data.MoleculeList.append(molecule)
    return


def readMeaningfulLine(file):
    readMore = True
    while (readMore):
        line = file.readline()
        index = line.find('!')
        line = line[:index]
        if (len(line.split()) != 0):
            return line


def readInertia(file):
    lines = file.readlines()
    numAtoms = int(lines[0])
    numRotors = int(lines[numAtoms + 1])
    pivots = []
    rotAtoms = []
    for i in range(numRotors):
        tokens = lines[numAtoms + 2 + 3 * i].split()
        pivots.append(int(tokens[0]))
        pivots.append(int(tokens[1]))
        atoms1 = []
        atoms2 = []
        atom1tokens = lines[numAtoms + 2 + 3 * i + 1].split()
        atom2tokens = lines[numAtoms + 2 + 3 * i + 2].split()
        for j in range(int(tokens[2])):
            atoms1.append(int(atom1tokens[j]))
        for j in range(int(tokens[3])):
            atoms2.append(int(atom2tokens[j]))
        rotAtoms.append(atoms1)
        rotAtoms.append(atoms2)
    return pivots, rotAtoms, numRotors


def readGeneralInertia(file, Mass):
    lines = file.readlines()
    rotors = []

    prevLevel = 1
    i = 0
    parentLevel = 0  # gives the index of parent in rotors

    for line in lines:
        tokens = line.split()
        if (len(tokens) == 0):
            continue
        level = int(tokens[0][1])

        if (level != 1):

            if (rotors[parentLevel].level == level - 2):
                parentLevel = i - 1
            elif (rotors[parentLevel].level >= level):
                jumps = rotors[parentLevel].level - level + 1
                for k in range(jumps):
                    parentLevel = rotors[parentLevel].parent
        symm = int(tokens[1])
        pivot2 = int(tokens[2])
        atomsList = []
        for atoms in tokens[3:]:
            atomsList.append(int(atoms))
        rotor = Rotor(atomsList, pivot2, level, symm, Mass)
        rotor.parent = parentLevel
        # print rotor.symm, rotor.pivot2, rotor.pivotAtom, rotor.atomsList
        rotors.append(rotor)
        i = i + 1
    return rotors


def readGeom(file):
    lines = file.readlines()
    lines.reverse()
# read geometries

    i = lines.index(
        "                          Input orientation:                          \n")
    geomLines = []
    stillRead = True
    k = i - 5
    while (stillRead):
        geomLines.append(lines[k])
        k = k - 1
        if (lines[k].startswith(" ------------")):
            stillRead = False

    geom = matrix(array(zeros((len(geomLines), 3), dtype=float)))
    Mass = matrix(array(zeros((len(geomLines), 1), dtype=float)))

    for j in range(len(geomLines)):
        line = geomLines[j]
        tokens = line.split()
        geom[j, 0] = double(tokens[3])
        geom[j, 1] = double(tokens[4])
        geom[j, 2] = double(tokens[5])
        if (int(tokens[1]) == 6):
            Mass[j] = 12.0
        if (int(tokens[1]) == 8):
            Mass[j] = 15.99491
        if (int(tokens[1]) == 1):
            Mass[j] = 1.00783
        if (int(tokens[1]) == 7):
            Mass[j] = 14.0031
        if (int(tokens[1]) == 17):
            Mass[j] = 34.96885
        if (int(tokens[1]) == 16):
            Mass[j] = 31.97207
        if (int(tokens[1]) == 9):
            Mass[j] = 18.99840
    return geom, Mass


def readGeomFc(file):
    lines = file.readlines()
    lines.reverse()
# read geometries

    i = lines.index(
        "                          Input orientation:                          \n")
    geomLines = []
    stillRead = True
    k = i - 5
    while (stillRead):
        geomLines.append(lines[k])
        k = k - 1
        if (lines[k].startswith(" ------------")):
            stillRead = False

    geom = matrix(array(zeros((len(geomLines), 3), dtype=float)))
    Mass = matrix(array(zeros((len(geomLines), 1), dtype=float)))

    for j in range(len(geomLines)):
        line = geomLines[j]
        tokens = line.split()
        geom[j, 0] = double(tokens[3])
        geom[j, 1] = double(tokens[4])
        geom[j, 2] = double(tokens[5])
        if (int(tokens[1]) == 6):
            Mass[j] = 12.0
        if (int(tokens[1]) == 8):
            Mass[j] = 15.99491
        if (int(tokens[1]) == 1):
            Mass[j] = 1.00783
        if (int(tokens[1]) == 7):
            Mass[j] = 14.0031
        if (int(tokens[1]) == 17):
            Mass[j] = 34.96885
        if (int(tokens[1]) == 16):
            Mass[j] = 31.97207
        if (int(tokens[1]) == 9):
            Mass[j] = 18.99840

# read force constants
    Fc = matrix(
        array(zeros((len(geomLines) * 3, len(geomLines) * 3), dtype=float)))
    i = lines.index(" Force constants in Cartesian coordinates: \n")
    fclines = []
    stillRead = True
    k = i - 1
    while (stillRead):
        fclines.append(lines[k])
        k = k - 1
        if (lines[k].startswith(" Force constants in internal coordinates:")):
            stillRead = False

    numRepeats = len(geomLines) * 3 / 5 + 1
    j = 0
    while j in range(numRepeats):
        i = 5 * j
        while i in range(5 * j, len(geomLines) * 3):
            line = fclines[(j * (len(geomLines) * 3 + 1) - j *
                            (j - 1) * 5 / 2) + i + 1 - 5 * j]
            tokens = line.split()
            k = 0
            while k in range(0, min(i + 1 - 5 * j, 5)):
                Fc[i, 5 * j + k] = float(tokens[k + 1].replace('D', 'E'))
                k = k + 1
            i = i + 1
        j = j + 1

    return geom, Mass, Fc


def readFc(file):
    lines = file.readlines()
    lines.reverse()
# read geometries

    i = lines.index(
        "                          Input orientation:                          \n")
    geomLines = []
    stillRead = True
    k = i - 5
    while (stillRead):
        geomLines.append(lines[k])
        k = k - 1
        if (lines[k].startswith(" ------------")):
            stillRead = False

    geom = matrix(array(zeros((len(geomLines), 3), dtype=float)))
    Mass = matrix(array(zeros((len(geomLines), 1), dtype=float)))

    for j in range(len(geomLines)):
        line = geomLines[j]
        tokens = line.split()
        geom[j, 0] = double(tokens[3])
        geom[j, 1] = double(tokens[4])
        geom[j, 2] = double(tokens[5])
        if (int(tokens[1]) == 6):
            Mass[j] = 12.0
        if (int(tokens[1]) == 8):
            Mass[j] = 15.99491
        if (int(tokens[1]) == 1):
            Mass[j] = 1.00783
        if (int(tokens[1]) == 7):
            Mass[j] = 14.0031
        if (int(tokens[1]) == 17):
            Mass[j] = 34.96885
        if (int(tokens[1]) == 16):
            Mass[j] = 31.97207
        if (int(tokens[1]) == 9):
            Mass[j] = 18.99840

# read force constants
    Fc = matrix(
        array(zeros((len(geomLines) * 3, len(geomLines) * 3), dtype=float)))
    i = lines.index(" Force constants in Cartesian coordinates: \n")
    fclines = []
    stillRead = True
    k = i - 1
    while (stillRead):
        fclines.append(lines[k])
        k = k - 1
        if (lines[k].startswith(" Force constants in internal coordinates:")):
            stillRead = False

    numRepeats = int(len(geomLines) * 3 / 5 + 1)
    j = 0
    while j in range(numRepeats):
        i = 5 * j
        while i in range(5 * j, len(geomLines) * 3):
            line = fclines[int((j * (len(geomLines) * 3 + 1) - j *
                            (j - 1) * 5 / 2) + i + 1 - 5 * j)]
            tokens = line.split()
            k = 0
            while k in range(0, min(i + 1 - 5 * j, 5)):
                Fc[i, 5 * j + k] = float(tokens[k + 1].replace('D', 'E'))
                k = k + 1
            i = i + 1
        j = j + 1

    return Fc


def read_freq(freq_file):
    with open(freq_file, 'r') as f:
        freq = []
        for line in f:
            # print(line)
            if re.search('Frequencies --', line):
                i = 0
                for token in line.split():
                    if i > 1: freq.append(float(token))
                    i += 1
        return freq


def printNormalModes(l, v, num, Mass):
    for i in range(num):
        for k in range(3):
            print('%18.3f ' % (sqrt(l[i * 3 + k]) * 337.0 / 6.5463e-02)),
        print

        for j in range(Mass.size):
            for m in range(3):
                for k in range(3):
                    print('%6.2f' % (v[j * 3 + k, i * 3 + m]))
                print(' ')
            print
        print


def readEnergy(file, string):
    com = file.read()
    if string == 'cbsqb3':
        Energy = re.search('CBS-QB3 \(0 K\)= ' +
                           ' \s*([\-0-9.]+)', com).group(1)
    elif string == 'g3':
        Energy = re.search('G3\(0 K\)= ' + ' \s*([\-0-9.]+)', com).group(1)

    elif string == 'ub3lyp':
        Energy = re.findall("E\(UB3LYP\) = " + "\s*([\-0-9.]+)", com)[-1]
    return float(Energy)


def readHFEnergy(fileName):
    file = open(fileName, 'r')
    Result = ""

    line = file.readline()

    startReading = False
    while line != "":

        if (line[:9] == r" 1\1\GINC"):
            startReading = True
        if (line[-4:-1] == r"\\@"):
            Result = Result + line[1:-1]
            break
        if (startReading):
            Result = Result + line[1:-1]
            # pdb.set_trace()
        line = file.readline()

    hf_5 = getNum(Result, r"\HF=")
    file.close()
    return hf_5


def getNum(Result, id):
    if (len(Result) <= 5):
        return 0
    fields = Result.split(id)
    resStr = fields[1]
    numstr = ""
    for i in range(len(resStr)):
        if resStr[i] == "\\":
            break
        numstr = numstr + resStr[i]
    return float(numstr)


def getAtoms(Mass):
    atoms = len(Mass) * ['']
    j = 0
    for m in Mass:
        if (int(m) == 12):
            atoms[j] = 'C'
        if (int(m) == 15):
            atoms[j] = 'O'
        if (int(m) == 1):
            atoms[j] = 'H'
        if (int(m) == 14):
            atoms[j] = 'N'
        if (int(m) == 34):
            atoms[j] = 'Cl'
        if (int(m) == 31):
            atoms[j] = 'S'
        if (int(m) == 18):
            atoms[j] = 'F'
        j = j + 1
    return atoms
