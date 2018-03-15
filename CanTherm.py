#!/usr/bin/env python
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

import sys
import readGeomFc
import pdb
import math
from numpy import *
from scipy import *
from constants import *
import kinetics



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
    out_file = 'output'
    rx = False

    # CBSQB3 E for H, N, O, C, P, Cl
    atomEcbsqb3 = {'H': -0.499818, 'N': -54.520543, 'O': -74.987624,
                   'C': -37.785385, 'P': -340.817186, 'Cl':-459.683658}
    # G3 E for H, N, O, C, P
    atomEg3 = {'H': -0.5010030, 'N': -54.564343, 'O': -75.030991, 'C': -37.827717,
               'P': -341.116432}
    # CCSD(T)-F12/cc-pVDZ-F12
    atomEccsdtf12 = {'H': -0.499811124128, 'N': -54.526406291655,
        'O': -74.995458316117, 'C': -37.788203485235, 'S': -397.663040369707}

    # CCSD(T)-F12/cc-pVTZ-F12
    atomEccsdt_f12_tz = {'H':-0.499946213243, 'N':-54.53000909621, 'O':-75.004127673424, 'C':-37.789862146471, 'S':-397.675447487865}

    # UB3LYP/cc-pVDZ
    atomEub3lyp = {'H': -0.501257936920, 'N': -54.4835759575, 'O':-75.0684969223,
        'C':-37.8519749084, 'S':-398.062555689 }

    # Experimental \DeltaH from NIST
    # https://webbook.nist.gov/cgi/cbook.cgi?Source=1984COX%2FWAG1B&Mask=1
    atomH = {'H': 52.103308, 'C': 171.29082008, 'N': 112.97335608,
        'O': 59.55551508, 'S': 66.24529302, 'Cl': 28.991666806}

    # Contains H + TC + SOC (spin orbital correction)
    #atomH = {'H': 50.62, 'N': 111.49, 'O': 58.163, 'C': 169.8147}

    # BAC for C-H C-C C=C C.TB.C  O-H  C-O C=O
    bondC = [-0.11, -0.3, -0.08, -0.64, 0.02, 0.33, 0.55]

    def __init__(self, input_filename=False, verbose=1):
        self.input_filename = input_filename
        self.verbose = verbose


    def run(self):
        inputFile = open(self.input_filename, 'r')

        # Read input file and calculate thermo
        readGeomFc.readInputFile(inputFile, self, verbose=self.verbose)

        # TODO Test this block and make sure see if it's even executing
        molecule = self.MoleculeList[0]
        H = molecule.Energy
        atoms = readGeomFc.getAtoms(molecule.Mass)
        atomsH = 0.0
        if molecule.Etype == 'cbsqb3':
            atomE = self.atomEcbsqb3
        if molecule.Etype == 'g3':
            atomE = self.atomEg3
        if molecule.Etype == 'ccsdtf12':
            atomE = self.atomEccsdtf12
        if molecule.Etype == 'DF-LUCCSD(T)-F12':
            atomE = self.atomEccsdt_f12_tz
        if molecule.Etype == 'ub3lyp':
            atomE = self.atomEub3lyp
        for atom in atoms:
            H -= atomE[atom]
            atomsH += self.atomH[atom]
        H = H * ha_to_kcal + atomsH
        #
        #     # if molecule.Etype == 'cbsqb3':
        #     #     b = 0
        #     #     for bonds in molecule.bonds:
        #     #         H += bonds * data.bondC[b]
        #     #         b += 1
        #     #
        #     # # TODO What's going on here
        #     # H += Thermal[i * len(Temp) + 0]
        #     # print('%12.2f' % H + '%12.2f' % Entropy[i * len(Temp) + 0])
        #     # for c in range(1, 8):
        #     #     print('%12.2f' % Cp[i * len(Temp) + c]),
        #     # print('\n')
        #

        # Calculate Kinetics
        if len(self.MoleculeList) == 1:
            self.write_output()
            return

        self.rx = kinetics.Reaction(self.MoleculeList[0], self.MoleculeList[1],
            self.Temp, tunneling="Wigner")
        self.rx.fit_arrhenius()
        # self.rx.print_arrhenius()

        # Write to output file
        self.write_output()




    def write_output(self):
        '''
        Write the output file to self.output
        '''

        self.oFile = open(self.out_file, 'w')
        oFile = self.oFile
        oFile.write(LICENSE + '\n\n\n\n')
        Temp = self.Temp

        # Write Molecular Data
        hr = '*'*80 + '\n\n'
        spacing = ' '*30
        oFile.write(hr + spacing + 'Molecule Properties\n\n' + hr)

        for molecule in self.MoleculeList:
            mol_header = 'Molecule %i\n' % (self.MoleculeList.index(molecule) + 1)
            mol_header += '-'*(len(mol_header)-1) + '\n\n'
            oFile.write(mol_header)
            molecule.print_properties(oFile)

        oFile.write('\n\n')

        # Write themo properties
        hr = '*'*80 + '\n\n'
        spacing = ' '*28
        oFile.write(hr + spacing + 'Thermodynamic Properties\n\n' + hr)

        for molecule in self.MoleculeList:
            mol_header = 'Molecule %i\n' % (self.MoleculeList.index(molecule) + 1)
            mol_header += '-'*(len(mol_header)-1) + '\n\n'
            oFile.write(mol_header)
            molecule.calculate_all_thermo(Temp,oFile)


        # Write Kinetics (if applicable)
        if self.rx != False:
            hr = '*'*80 + '\n\n'
            spacing = ' '*31
            oFile.write(hr + spacing + 'Kinetic Properties\n\n' + hr)


            self.rx.print_properties(oFile)

        # Print goodbye
        oFile.write('\n\n\n')
        goodbye = ' '*19 + 'CANTHERM OUTPUT COMPLETE. HAVE A GOOD DAY\n'
        goodbye = '*'*80 + '\n\n' + goodbye + '\n'  + '*'*80 + '\n'
        oFile.write(goodbye)
        oFile.close()


if __name__ == "__main__":
    # main()
    data = CanTherm(sys.argv[1])
    data.run()
