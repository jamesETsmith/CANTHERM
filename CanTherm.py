#!/usr/bin/env python

import sys
import readGeomFc
import pdb
import math
from numpy import *
from scipy import *
from constants import *


def wigner_correction(T, freq, scale):
    # see doi:10.1103/PhysRev.40.749 and doi:10.1039/TF9595500001
    return (1.0 + 1.0 / 24.0 * (h * abs(freq) * scale * c_in_cm / (T * kb) )**2)

class CanTherm:
    CalcType = ''
    ReacType = ''
    Temp = []
    MoleculeList = []
    Entropy = []
    Thermal = []
    Cp = []
    Parition = []
    scale = 0.0
    # CBSQB3 E for H, N, O, C, P
    atomEcbsqb3 = {'H': -0.499818, 'N': -54.520543, 'O': -74.987624,
                   'C': -37.785385, 'P': -340.817186}
    # G3 E for H, N, O, C, P
    atomEg3 = {'H': -0.5010030, 'N': -54.564343, 'O': -75.030991, 'C': -37.827717,
               'P': -341.116432}
    # CCSD(T)-F12/cc-pVDZ-F12
    atomEccsdtf12 = {'H': -0.499811124128, 'N': -54.526406291655,
        'O': -74.995458316117, 'C': -37.788203485235, 'S': -397.663040369707}

    # UB3LYP/cc-pVDZ
    atomEub3lyp = {'H': -0.501257936920, 'N': -54.4835759575, 'O':-75.0684969223,
        'C':-37.8519749084, 'S':-398.062555689 }

    # expt H contains H + TC + SOC (spin orbital correction)
    atomH = {'H': 50.62, 'N': 111.49, 'O': 58.163, 'C': 169.8147}
    # BAC for C-H C-C C=C C.TB.C  O-H  C-O C=O
    bondC = [-0.11, -0.3, -0.08, -0.64, 0.02, 0.33, 0.55]


def main():
    data = CanTherm()
    inputFile = open(sys.argv[1], 'r')
    oFile = open('output', 'w')
    readGeomFc.readInputFile(inputFile, data)

    data.Entropy = len(data.MoleculeList) * len(data.Temp) * [0.0]
    data.Cp = len(data.MoleculeList) * len(data.Temp) * [0.0]
    data.Thermal = len(data.MoleculeList) * len(data.Temp) * [0.0]
    data.Partition = len(data.MoleculeList) * len(data.Temp) * [1.0]
    Entropy = data.Entropy
    Cp = data.Cp
    Thermal = data.Thermal
    Partition = data.Partition

    for i in range(len(data.MoleculeList)):
        molecule = data.MoleculeList[i]
        oFile.write('Molecule ' + str(i + 1) + ':\n')
        oFile.write('-----------\n\n')
        molecule.printData(oFile)

        oFile.write('\nThermodynamic Data\n')

        Temp = data.Temp
        # translation
        (ent, cp, dh) = molecule.getTranslationThermo(oFile, data.Temp)
        for j in range(len(Temp)):
            Entropy[i * len(Temp) + j] = Entropy[i * len(Temp) + j] + ent[j]
            Cp[i * len(Temp) + j] = Cp[i * len(Temp) + j] + cp[j]

        # vibrational
        (ent, cp, dh, q) = molecule.getVibrationalThermo(
            oFile, data.Temp, data.scale)
        for j in range(len(Temp)):
            Entropy[i * len(Temp) + j] = Entropy[i * len(Temp) + j] + ent[j]
            Cp[i * len(Temp) + j] = Cp[i * len(Temp) + j] + cp[j]
            Thermal[i * len(Temp) + j] = Thermal[i * len(Temp) + j] + dh[j]

        # Internal rotational
        if molecule.numRotors != 0:
            (ent, cp, dh, q) = molecule.getIntRotationalThermo_Q(oFile, data.Temp)
            for j in range(len(Temp)):
                Entropy[i * len(Temp) + j] = Entropy[i * len(Temp) + j] + ent[j]
                Cp[i * len(Temp) + j] = Cp[i * len(Temp) + j] + cp[j]
                Thermal[i * len(Temp) + j] = Thermal[i * len(Temp) + j] + dh[j]

        # External rotational
        (ent, cp, dh) = molecule.getExtRotationalThermo(oFile, data.Temp)
        for j in range(len(Temp)):
            Entropy[i * len(Temp) + j] = Entropy[i * len(Temp) + j] + ent[j]
            Cp[i * len(Temp) + j] = Cp[i * len(Temp) + j] + cp[j]
            Thermal[i * len(Temp) + j] = Thermal[i * len(Temp) + j] + dh[j]

        for j in range(len(Temp)):
            Entropy[i * len(Temp) + j] = Entropy[i * len(Temp) +
                                                 j] + R_kcal * math.log(molecule.nelec)

        H = molecule.Energy
        atoms = readGeomFc.getAtoms(molecule.Mass)
        atomsH = 0.0
        if molecule.Etype == 'cbsqb3':
            atomE = data.atomEcbsqb3
        if molecule.Etype == 'g3':
            atomE = data.atomEg3
        if molecule.Etype == 'ccsdtf12':
            atomE = data.atomEccsdtf12
        if molecule.Etype == 'ub3lyp':
            atomE = data.atomEub3lyp
        for atom in atoms:
            H -= atomE[atom]
            atomsH += data.atomH[atom]
        H = H * ha_to_kcal + atomsH

        if molecule.Etype == 'cbsqb3':
            b = 0
            for bonds in molecule.bonds:
                H += bonds * data.bondC[b]
                b += 1

        H += Thermal[i * len(Temp) + 0]
        print('%12.2f' % H + '%12.2f' % Entropy[i * len(Temp) + 0])
        for c in range(1, 8):
            print('%12.2f' % Cp[i * len(Temp) + c]),
        print('\n')

    if len(data.MoleculeList) == 1:
        return

    # fit the rate coefficient
    A = matrix(zeros((len(Temp), 3), dtype=float))
    y = matrix(zeros((len(Temp), 1), dtype=float))

    rate = [0.0] * len(Temp)
    kappa = []
    for j in range(len(Temp)):
        if (data.ReacType == 'Unimol'):
            rate[j] = (kb * Temp[j] / h) * math.exp((Entropy[len(Temp) + j] - Entropy[j]) / R) * \
                math.exp(-(data.MoleculeList[1].Energy - data.MoleculeList[0].Energy)
                         * ha_to_kcal * 1.0e3 / R_kcal / Temp[j])

            kappa.append( wigner_correction(Temp[j], data.MoleculeList[1].imagFreq, data.scale ) )

            rate[j] *= kappa[j]
            A[j, :] = mat([1.0, math.log(Temp[j]), -1.0 / R_kcal / Temp[j]])
            y[j] = log(rate[j])

    b = linalg.inv(transpose(A) * A) * (transpose(A) * y)

    # Write the reaction rate data
    oFile.write('\n\nRate Data\n')
    oFile.write('r = A*(T/1000)^n*exp(-Ea/R/T)' + '%12.2e' % (exp(b[0]) * 1000.0**float(b[1])) + '%10.2f' % b[1] + '%12.2f' % (b[2] / 1.0e3) +
                '\n') #TODO what is this?
    oFile.write('r = A*T^n*exp(-Ea/R/T)' + '%12.2e' %
                (exp(b[0])) + '%10.2f' % b[1] + '%12.2f' % (b[2] / 1.0e3) + '\n\n')

    oFile.write('%12s' % 'Temp. (K)' + '%16s' % 'Rate (s^-1)'+ '%16s' % 'FitRate (s^-1)'  + '%12s\n' % 'Kappa')

    for j in range(len(Temp)):
        fitrate = exp(b[0]) * Temp[j]**float(b[1]) * \
            exp(-b[2] / R_kcal / Temp[j])
        oFile.write('%12.2f' % Temp[j] + '%16.2e' %
                    rate[j] + '%16.2e' % fitrate + '%12.4f\n' % kappa[j])
    oFile.write('\n\n')
    oFile.close()


if __name__ == "__main__":
    main()
