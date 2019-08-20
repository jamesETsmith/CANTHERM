#!/usr/bin/env python
# takes the file as an argument and then reads the geometry and lower triangular part
# of the force constant matrix and the Atomic mass vector
'''
********************************************************************************
*                                                                              *
*    CANTHERM                                                                  *
*                                                                              *
*    Copyright (C) 2018, Sandeep Sharma and James E. T. Smith                  *
*                                                                              *
*    This program is free software: you can redistribute it and/or modify      *
*    it under the terms of the GNU General Public License as published by      *
*    the Free Software Foundation, either version 3 of the License, or         *
*    (at your option) any later version.                                       *
*                                                                              *
*    This program is distributed in the hope that it will be useful,           *
*    but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*    GNU General Public License for more details.                              *
*                                                                              *
*    You should have received a copy of the GNU General Public License         *
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
*                                                                              *
********************************************************************************
'''

import os, re
import numpy as np

from cantherm import rotor
from cantherm import molecule


def readInputFile(file, data, verbose):

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
            print(
                'ReactionType is not specified in the. Either Unimol or Bimol should be specified after keyword Read'
            )
            exit()
        line = readMeaningfulLine(file)


# read Temperature range
    if line.split()[0].upper() == 'TLIST':
        line = readMeaningfulLine(file)
        numTemp = int(line.split()[0])
        data.Temp = []
        i = 0
        while i < numTemp:
            line = readMeaningfulLine(file)
            tokens = line.split()
            i = i + len(tokens)
            for j in tokens:
                data.Temp.append(float(j))
        if len(data.Temp) > numTemp:
            print(len(data.Temp), numTemp)
            print('More Temperatures than ', numTemp, ' are specified')

    elif line.split()[0].upper() == 'TRANGE':
        line = readMeaningfulLine(file)
        tokens = line.split()
        T0 = float(tokens[0])
        dT = float(tokens[1])
        numTemp = int(tokens[2])
        data.Temp = []
        for i in range(numTemp):
            data.Temp.append(T0 + i * dT)
    else:
        print(
            'Temperaure information not given either use keyword Tlist or Trange'
        )
        exit()

    # read scaling factor
    line = readMeaningfulLine(file)
    if (line.split()[0].upper() != 'SCALE'):
        print('Give a scaling factor')
        exit(0)
    data.scale = float(line.split()[1])

    # Read root directory (if included)
    line = readMeaningfulLine(file)
    if (line.split()[0].upper() == 'ROOTDIR'):
        data.dir = line.split()[1] + '/'
    else:
        print('Specify root directory')
        exit(0)

    # read molecular data
    if (data.CalcType == 'Thermo'):
        numMol = 1
    elif data.CalcType == 'Reac' and data.ReacType == 'Unimol':
        numMol = 2
    elif data.CalcType == 'Reac' and data.ReacType == 'Bimol':
        numMol = 3

    for i in range(numMol):
        if data.ReacType == 'Unimol' and i == 1 or data.ReacType == 'Bimol' and i == 2:
            molec = molecule.Molecule(file, True, data.scale, verbose,
                                      data.dir)
        else:
            molec = molecule.Molecule(file, False, data.scale, verbose,
                                      data.dir)
        data.MoleculeList.append(molec)
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
        rotor = rotor.Rotor(atomsList, pivot2, level, symm, Mass)
        rotor.parent = parentLevel
        # print rotor.symm, rotor.pivot2, rotor.pivotAtom, rotor.atomsList
        rotors.append(rotor)
        i = i + 1
    return rotors


def readGeom(file):
    # print("HEREHEREHERE") # TODO
    lines = file.readlines()
    lines.reverse()
    # read geometries

    i = lines.index(
        "                          Input orientation:                          \n"
    )
    geomLines = []
    stillRead = True
    k = i - 5
    while (stillRead):
        geomLines.append(lines[k])
        k = k - 1
        if (lines[k].startswith(" ------------")):
            stillRead = False

    geom = np.zeros((len(geomLines), 3))
    Mass = np.zeros((len(geomLines)))

    for j in range(len(geomLines)):
        line = geomLines[j]
        tokens = line.split()
        geom[j, 0] = float(tokens[3])
        geom[j, 1] = float(tokens[4])
        geom[j, 2] = float(tokens[5])
        if (int(tokens[1]) == 6):
            Mass[j] = 12.0
        elif (int(tokens[1]) == 8):
            Mass[j] = 15.99491
        elif (int(tokens[1]) == 1):
            Mass[j] = 1.00783
        elif (int(tokens[1]) == 7):
            Mass[j] = 14.0031
        elif (int(tokens[1]) == 17):
            Mass[j] = 34.96885
        elif (int(tokens[1]) == 16):
            Mass[j] = 31.97207
        elif (int(tokens[1]) == 9):
            Mass[j] = 18.99840

        elif (int(tokens[1]) == 13):
            Mass[j] = 26.98153853
        elif (int(tokens[1]) == 14):
            Mass[j] = 27.97692653465
        elif (int(tokens[1]) == 15):
            Mass[j] = 30.97376199842

        elif (int(tokens[1]) == 18):
            Mass[j] = 35.96754510
        elif (int(tokens[1]) == 19):
            Mass[j] = 38.9637064864
        elif (int(tokens[1]) == 20):
            Mass[j] = 39.962590863
        elif (int(tokens[1]) == 21):
            Mass[j] = 44.95590828
        elif (int(tokens[1]) == 22):
            Mass[j] = 45.95262772
        elif (int(tokens[1]) == 23):
            Mass[j] = 49.94715601
        elif (int(tokens[1]) == 24):
            Mass[j] = 49.94604183
        elif (int(tokens[1]) == 25):
            Mass[j] = 54.93804391
        elif (int(tokens[1]) == 26):
            Mass[j] = 53.93960899
        elif (int(tokens[1]) == 27):
            Mass[j] = 58.93319429
        elif (int(tokens[1]) == 28):
            Mass[j] = 57.93534241
        elif (int(tokens[1]) == 29):
            Mass[j] = 62.92959772
        elif (int(tokens[1]) == 30):
            Mass[j] = 63.92914201
        elif (int(tokens[1]) == 31):
            Mass[j] = 68.9255735
        elif (int(tokens[1]) == 32):
            Mass[j] = 69.92424875
        elif (int(tokens[1]) == 33):
            Mass[j] = 74.92159457
        elif (int(tokens[1]) == 34):
            Mass[j] = 73.922475934
        elif (int(tokens[1]) == 35):
            Mass[j] = 78.9183376
        elif (int(tokens[1]) == 36):
            Mass[j] = 77.9203649
        elif (int(tokens[1]) == 37):
            Mass[j] = 84.9117897379
        elif (int(tokens[1]) == 38):
            Mass[j] = 83.9134191
        elif (int(tokens[1]) == 39):
            Mass[j] = 88.9058403
        elif (int(tokens[1]) == 40):
            Mass[j] = 89.9046977
        elif (int(tokens[1]) == 41):
            Mass[j] = 92.906373
        elif (int(tokens[1]) == 42):
            Mass[j] = 91.90680796
        elif (int(tokens[1]) == 43):
            Mass[j] = 96.9063667
        elif (int(tokens[1]) == 44):
            Mass[j] = 95.90759025
        elif (int(tokens[1]) == 45):
            Mass[j] = 102.905498
        elif (int(tokens[1]) == 46):
            Mass[j] = 101.9056022
        elif (int(tokens[1]) == 47):
            Mass[j] = 106.9050916
        elif (int(tokens[1]) == 48):
            Mass[j] = 105.9064599
        elif (int(tokens[1]) == 49):
            Mass[j] = 112.90406184
        elif (int(tokens[1]) == 50):
            Mass[j] = 111.90482387
        elif (int(tokens[1]) == 51):
            Mass[j] = 120.903812
        elif (int(tokens[1]) == 52):
            Mass[j] = 119.9040593
        elif (int(tokens[1]) == 53):
            Mass[j] = 126.9044719
        elif (int(tokens[1]) == 54):
            Mass[j] = 123.905892
        elif (int(tokens[1]) == 55):
            Mass[j] = 132.905451961
        elif (int(tokens[1]) == 56):
            Mass[j] = 129.9063207
        elif (int(tokens[1]) == 57):
            Mass[j] = 137.9071149
        elif (int(tokens[1]) == 58):
            Mass[j] = 135.90712921
        elif (int(tokens[1]) == 59):
            Mass[j] = 140.9076576
        elif (int(tokens[1]) == 60):
            Mass[j] = 141.907729
        elif (int(tokens[1]) == 61):
            Mass[j] = 144.9127559
        elif (int(tokens[1]) == 62):
            Mass[j] = 143.9120065
        elif (int(tokens[1]) == 63):
            Mass[j] = 150.9198578
        elif (int(tokens[1]) == 64):
            Mass[j] = 151.9197995
        elif (int(tokens[1]) == 65):
            Mass[j] = 158.9253547
        elif (int(tokens[1]) == 66):
            Mass[j] = 155.9242847
        elif (int(tokens[1]) == 67):
            Mass[j] = 164.9303288
        elif (int(tokens[1]) == 68):
            Mass[j] = 161.9287884
        elif (int(tokens[1]) == 69):
            Mass[j] = 168.9342179
        elif (int(tokens[1]) == 70):
            Mass[j] = 167.9338896
        elif (int(tokens[1]) == 71):
            Mass[j] = 174.9407752
        elif (int(tokens[1]) == 72):
            Mass[j] = 173.9400461
        elif (int(tokens[1]) == 73):
            Mass[j] = 179.9474648
        elif (int(tokens[1]) == 74):
            Mass[j] = 179.9467108
        elif (int(tokens[1]) == 75):
            Mass[j] = 184.9529545
        elif (int(tokens[1]) == 76):
            Mass[j] = 183.9524885
        elif (int(tokens[1]) == 77):
            Mass[j] = 190.9605893
        elif (int(tokens[1]) == 78):
            Mass[j] = 189.9599297
        elif (int(tokens[1]) == 79):
            Mass[j] = 196.96656879
        elif (int(tokens[1]) == 80):
            Mass[j] = 195.9658326
        elif (int(tokens[1]) == 81):
            Mass[j] = 202.9723446
        elif (int(tokens[1]) == 82):
            Mass[j] = 203.973044
        elif (int(tokens[1]) == 83):
            Mass[j] = 208.9803991
        elif (int(tokens[1]) == 84):
            Mass[j] = 208.9824308
        elif (int(tokens[1]) == 85):
            Mass[j] = 209.9871479
        elif (int(tokens[1]) == 86):
            Mass[j] = 210.9906011
        else:
            print("ATOM NOT FOUND WITH TOKEN: %s" % tokens[1])

    # Read the bonds
    bond_ind = lines.index(
        '                           !    Initial Parameters    !\n')
    bond_ind -= 5  # Skip heading
    reading_bonds = True
    bonds = []
    while reading_bonds:
        if re.search(' ! R\d', lines[bond_ind]):
            s = re.search('\(([^)]+)', lines[bond_ind]).group(1).split(',')
            s = [int(s[0]) - 1, int(s[1]) - 1]
            bonds.append(s)
            # print(s) # TODO
        if re.search(' ! A\d', lines[bond_ind]):
            reading_bonds = False
            break
        if re.search(' ! D\d', lines[bond_ind]):
            reading_bonds = False
            break
        bond_ind -= 1
    # exit(0) # TODO
    return geom, Mass, bonds


def readGeomFc(file):
    lines = file.readlines()
    lines.reverse()
    # read geometries

    i = lines.index(
        "                          Input orientation:                          \n"
    )
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
        geom[j, 0] = float(tokens[3])
        geom[j, 1] = float(tokens[4])
        geom[j, 2] = float(tokens[5])
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
            line = fclines[(j * (len(geomLines) * 3 + 1) - j * (j - 1) * 5 / 2)
                           + i + 1 - 5 * j]
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
        "                          Input orientation:                          \n"
    )
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
        geom[j, 0] = float(tokens[3])
        geom[j, 1] = float(tokens[4])
        geom[j, 2] = float(tokens[5])
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
    tokens = file.split('.')
    efile = open(file, 'r')
    com = efile.read()

    # Gaussian File
    if tokens[-1] == 'log' or tokens[-1] == 'out':
        if string == 'cbsqb3':
            Energy = re.search('CBS-QB3 \(0 K\)= ' + ' \s*([\-0-9.]+)',
                               com).group(1)
        elif string == 'g3':
            Energy = re.search('G3\(0 K\)= ' + ' \s*([\-0-9.]+)', com).group(1)

        elif string == 'ub3lyp':
            Energy = re.findall("E\(UB3LYP\) = " + "\s*([\-0-9.]+)", com)[-1]

        elif string == 'RM062X':
            Energy = re.findall("E\(RM062X\) = " + "\s*([\-0-9.]+)", com)[-1]
        elif string == 'UCCSD(T)':
            Energy = re.findall(" E\(Corr\)= " + "\s*([\-0-9.]+)", com)[-1]
            print(Energy)
        else:
            raise ValueError("Gaussian method {} is not supported.".format(string))

    # Molpro File
    elif tokens[-1] == 'res':
        if string == 'DF-LUCCSD(T)-F12':
            Energy = re.search('DF-LUCCSD\(T\)-F12\/cc-pVTZ-F12 energy=' + \
                                 ' \s*([\-0-9.]+)', com).group(1)

        elif string == 'UCCSD(T)-F12':
            Energy = re.search('UCCSD\(T\)-F12\/cc-pVTZ-F12 energy=' + \
                                 ' \s*([\-0-9.]+)', com).group(1)

        elif string == 'RHF-LRMP2':
            Energy = re.search('RHF-LRMP2 STATE 1.1 Energy' + '\s*([\-0-9.]+)',
                               com).group(1)

        elif string == 'RHF-RMP2':
            Energy = re.search('RHF-RMP2 energy' + '\s*([\-0-9.]+)',
                               com).group(1)

        elif string == 'RHF-UCCSD-F12a':
            Energy = re.search('RHF-UCCSD-F12a energy' + '\s*([\-0-9.]+)',
                               com).group(1)

        elif string == 'LUCCSD-F12a':
            E_corr = re.search(
                'LUCCSD-F12a correlation energy' + '\s*([\-0-9.]+)',
                com).group(1)
            E_ref = re.search('New reference energy' + '\s*([\-0-9.]+)',
                              com).group(1)
            Energy = float(E_ref) + float(E_corr)
        else:
            raise ValueError("MOLPRO method {} is not supported.".format(string))

    efile.close()
    return float(Energy)


def readMOEnergies(filename):

    aocc = " Alpha  occ. eigenvalues --"
    avir = " Alpha virt. eigenvalues --"
    bocc = "  Beta  occ. eigenvalues --"
    bvir = "  Beta virt. eigenvalues --"
    end_str = " **********************************************************************"

    with open(filename, 'r') as f:
        lines = f.readlines()

        reading_mos = False
        raw_mos = []

        # Read MOs
        for l in range(1, len(lines)):
            i = -l
            # Check for strings to start reading MOs
            # print(lines[i])
            if not reading_mos:
                if bvir in lines[i]:
                    reading_mos = True
                elif bocc in lines[i]:
                    reading_mos = True
                elif avir in lines[i]:
                    reading_mos = True
                elif aocc in lines[i]:
                    reading_mos = True

                if reading_mos:
                    if len(lines[i]) > 1:
                        # print(lines[i])
                        for mo in lines[i].split()[5:]:
                            raw_mos.append(float(mo))

            else:
                # Test to stop and end reading of MOs
                if end_str in lines[i]:
                    reading_mos = False
                    break

                # Read data
                if len(lines[i]) > 1:
                    # print(len(lines[i]))
                    # print(lines[i])
                    for mo in lines[i].split()[5:]:
                        raw_mos.append(float(mo))

        mos = np.sort(np.array(raw_mos))
        homo = mos[mos < 0.0].max()
        lumo = mos[mos > 0.0].min()
        # print(mos, homo, lumo)
        # print(mos)
        return mos, homo, lumo


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
    atoms = Mass.size * ['']
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


def read_dipole(dipole_file):
    '''
    Returns the dipole and its magnitude in Debye.
    '''

    with open(dipole_file, 'r') as f:
        lines = f.readlines()
        n_lines = len(lines)

        magn = 0
        dip = np.zeros((3, ))

        for i in range(n_lines):
            if ' Electric dipole moment (dipole orientation)' in lines[
                    n_lines - i - 1]:
                magn = float(lines[n_lines - i + 2].split()[2].replace(
                    'D', 'E'))
                dip[0] = float(lines[n_lines - i + 3].split()[2].replace(
                    'D', 'E'))
                dip[1] = float(lines[n_lines - i + 4].split()[2].replace(
                    'D', 'E'))
                dip[2] = float(lines[n_lines - i + 5].split()[2].replace(
                    'D', 'E'))
                break

        # print(magn, dip) TODO
        return magn, dip  # in Debye
