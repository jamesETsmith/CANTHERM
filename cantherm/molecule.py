#!/usr/bin/env python
"""
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
"""

import os
import numpy as np

from cantherm import readGeomFc
from cantherm.constants import *
from cantherm import rotor

la = np.linalg


class Molecule:
    """
    Attributes:

        input_file (open file object): Filename for cantherm input

        TS (boolean): Is the molecule a transition state.

        imagFreq (float): Imaginary frequency of transition state in
            :math:`cm^{-1}`.

        scale (float): Frequency scale factor.

        verbose (int): Controls output. 0 means no output at all.

        Freq ([float]): Harmonic frequencies.

        hindFreq ([float]): Harmonic frequencies for internal rotors.

        numRotors (int): Number of hindered rotors.

        rotors ([Rotor]): List of :class:`Rotor` objects specified in the input.

        bonds ([[int]]): 2D list of the bonds (atoms are the atom number - 1).

        Energy (float): Energy from log file in hartrees.

        Etype (string): Type of energy calculation performed.

        geom (np.array): 2D numpy array of floats that holds the coords in
            angstroms.

        Mass (np.array): 1D numpy array of floats for the masses, index matches
            atom number.

        e_file (string): Path to the energy file (if one exists)

        dir (string): Root directory for all FILEs listed in input.

        name (string): Molecule label.

        calc_complete (bool): True if calculation completed correctly.

    """

    def __init__(self, in_file, isTS, scale, verbose, root_dir, txt_input=True):
        self.input_file = in_file
        self.TS = isTS
        self.scale = scale
        self.verbose = verbose
        self.dir = root_dir

        # Initialize other attributes
        self.Freq = []
        self.hindFreq = []
        self.rotors = []
        self.bonds = []
        self.Etype = ""
        self.label = ""
        self.calc_complete = True

        # Assign rest of attributes
        if txt_input:
            self.read_input()
        else:
            self.read_log(in_file)  # TODO this is not a great way to do this

    def read_input(self):
        """
        Reads the input file and associated log (or other output format) files.
        """

        file = self.input_file
        line = readGeomFc.readMeaningfulLine(file)
        self.label = line.split()[-1]

        # read linearity
        line = readGeomFc.readMeaningfulLine(file)
        self.linearlity = line.split()[0].upper

        # read geometry
        line = readGeomFc.readMeaningfulLine(file)
        tokens = line.split()
        if tokens[0].upper() != "GEOM":
            print(tokens)
            print("Geom keyword not found in the input file")
            exit()

        if len(tokens) == 1:
            # geometry is given following the GEOM keyword
            line = readGeomFc.readMeaningfulLine(file)
            numAtoms = int(line.split()[0])
            for i in range(numAtoms):
                line = readGeomFc.readMeaningfulLine(file)
                tokens = line.split()
                self.geom[j, 0] = double(tokens[3])
                self.geom[j, 1] = double(tokens[4])
                self.geom[j, 2] = double(tokens[5])
                if int(tokens[1]) == 6:
                    self.Mass[j] = 12.00000
                if int(tokens[1]) == 8:
                    self.Mass[j] = 15.99491
                if int(tokens[1]) == 1:
                    self.Mass[j] = 1.00783
                if int(tokens[1]) == 7:
                    self.Mass[j] = 14.0031
                if int(tokens[1]) == 17:
                    self.Mass[j] = 34.96885
                if int(tokens[1]) == 16:
                    self.Mass[j] = 31.97207
                if int(tokens[1]) == 9:
                    self.Mass[j] = 18.99840

        # read geometry from the file
        if tokens[1].upper() == "FILE":
            if self.verbose > 0:
                print("Reading Geometry from the file: ", self.dir + tokens[2])
            geomFile = open(self.dir + tokens[2], "r")
            (self.geom, self.Mass, self.bonds) = readGeomFc.readGeom(geomFile)
            geomFile.close()
            # print self.geom
        else:
            print(
                "Either give geometry or give keyword File followed by the file containing geometry data"
            )
            exit()

        self.calculateMomInertia()

        # read force constant or frequency data
        line = readGeomFc.readMeaningfulLine(file)
        tokens = line.split()
        if tokens[0].upper() == "FORCEC" and tokens[1].upper() == "FILE":
            fcfile = open(self.dir + tokens[2], "r")
            self.Fc = readGeomFc.readFc(fcfile)

            for i in range(0, 3 * self.Mass.size):
                for j in range(i, 3 * self.Mass.size):
                    self.Fc[i, j] = self.Fc[j, i]

        elif tokens[0].upper() == "FREQ" and len(tokens) == 1:
            line = readGeomFc.readMeaningfulLine(file)
            numFreq = int(line.split()[0])
            i = 0
            while i < numFreq:
                line = readGeomFc.readMeaningfulLine(file)
                tokens = line.split()
                i = i + len(tokens)
                for j in tokens:
                    if float(j) < 0 and self.TS:
                        print("Imaginary Freq", j)
                        self.imagFreq = float(j)
                        continue
                    self.Freq.append(float(j))

            if len(self.Freq) > numFreq:
                print("More frequencies than ", numFreq, " are specified")

        elif tokens[0].upper() == "FREQ" and tokens[1].upper() == "FILE":
            freq_raw = readGeomFc.read_freq(self.dir + tokens[2])
            for freq in freq_raw:
                if freq < 0 and self.TS:
                    self.imagFreq = freq
                else:
                    self.Freq.append(freq)

        else:
            print("Frequency information cannot be read, check input file again")
            exit()

        # read energy
        line = readGeomFc.readMeaningfulLine(file)
        tokens = line.split()
        if tokens[0].upper() != "ENERGY":
            print("Energy information not given")
            exit()
        if tokens[1].upper() == "FILE":
            if self.verbose > 0:
                print("Reading energy from file: ", self.dir + tokens[2])
            if tokens[3].upper() == "CBS-QB3":
                self.Etype = "cbsqb3"
            elif tokens[3].upper() == "G3":
                self.Etype = "g3"
            elif tokens[3].upper() == "UB3LYP":
                self.Etype = "ub3lyp"
            elif tokens[3].upper() == "RM062X":
                self.Etype = "RM062X"
            elif tokens[3].upper() == "DF-LUCCSD(T)-F12":
                self.Etype = "DF-LUCCSD(T)-F12"
            elif tokens[3].upper() == "UCCSD(T)-F12":
                self.Etype = "UCCSD(T)-F12"
            self.Energy = readGeomFc.readEnergy(self.dir + tokens[2], self.Etype)
            self.e_file = self.dir + tokens[2]
            if self.verbose > 0:
                print(self.Energy, self.Etype)
        elif len(tokens) == 3:
            self.e_file = None
            self.Energy = float(tokens[1])
            if tokens[3].upper() == "CBS-QB3":
                self.Etype = "cbsqb3"
            elif tokens[3].upper() == "G3":
                self.Etype = "g3"
            elif tokens[3].upper() == "UB3LYP":
                self.Etype = "ub3lyp"
            elif tokens[3].upper() == "RM062X":
                self.Etype = "RM062X"
            elif tokens[3].upper() == "DF-LUCCSD(T)-F12":
                self.Etype = "DF-LUCCSD(T)-F12"
            elif tokens[3].upper() == "UCCSD(T)-F12":
                self.Etype = "UCCSD(T)-F12"
            if self.verbose > 0:
                print(self.Etype.upper(), " Energy: ", self.Energy)
        else:
            print("Cannot read the Energy")
            exit()

        # read external symmetry
        line = readGeomFc.readMeaningfulLine(file)
        if line.split()[0].upper() != "EXTSYM":
            print("Extsym keyword required")
            exit()
        self.extSymm = float(line.split()[1])

        # read electronic degeneracy
        line = readGeomFc.readMeaningfulLine(file)
        if line.split()[0].upper() != "NELEC":
            print("Nelec keyword required")
            exit()
        self.nelec = int(line.split()[1])

        # read rotor information
        #         line = readGeomFc.readMeaningfulLine(file)
        #         if (line.split()[0].upper() != 'ROTORS'):
        #             print('Rotors keyword required')
        #             exit()
        #         self.numRotors = int(line.split()[1])
        #         if self.numRotors == 0:
        #             self.rotors = []
        #             if (len(self.Freq) == 0):
        #                 # calculate frequencies from force constant
        #                 self.getFreqFromFc()
        #             return
        #
        #         rotorFile = line.split()[2]
        #         inertiaFile = open(rotorFile, 'r')
        #         # print self.Mass
        #         (self.rotors) = readGeomFc.readGeneralInertia(inertiaFile, self.Mass)
        #         if len(self.rotors) - 1 != self.numRotors:
        #             print("The number of rotors specified in file, ",
        #                   rotorFile, ' is different than, ', self.numRotors)
        #
        #         if (len(self.Freq) == 0):
        #             # calculate frequencies from force constant
        #             self.getFreqFromFc()
        #
        #
        # # read potential information for rotors
        #         line = readGeomFc.readMeaningfulLine(file)
        #         tokens = line.split()
        #         if tokens[0].upper() != 'POTENTIAL':
        #             print('No information for potential given')
        #             exit()
        #
        #         if tokens[1].upper() == 'SEPARABLE':
        #             if tokens[2].upper() == 'FILES':
        #                 line = readGeomFc.readMeaningfulLine(file)
        #                 tokens = line.split()
        #                 if len(tokens) != self.numRotors:
        #                     print('give a separate potential file for each rotor')
        #                 for files in tokens:
        #                     Kcos = []
        #                     Ksin = []
        #                     harmonic = Harmonics(5, Kcos, Ksin)
        #                     harmonic.fitPotential(files)
        #                     self.Harmonics.append(harmonic)
        #
        #             elif tokens[2].upper() == 'HARMONIC':
        #                 for i in range(self.numRotors):
        #                     line = readGeomFc.readMeaningfulLine(file)
        #                     numFit = int(line.split()[0])
        #                     Ksin = []
        #                     Kcos = []
        #                     for i in range(numFit):
        #                         line = readGeomFc.readMeaningfulLine(file)
        #                         tokens = line.split()
        #                         Kcos.append(tokens[0])
        #                         Ksin.append(tokens[0])
        #                     harmonic = Harmonics(numFit, Kcos, Ksin)
        #                     self.Harmonics.append(harmonic)
        #
        #         elif tokens[1].upper() == 'NONSEPARABLE':
        #             line = readGeomFc.readMeaningfulLine(file)
        #             self.potentialFile = line.split()[0]

        # Rotors
        line = readGeomFc.readMeaningfulLine(file)
        tokens = line.split()

        if tokens[0].upper() != "ROTORS":
            print("No information for the potential given.\nExiting...")
            exit()

        self.numRotors = int(line.split()[1])

        if self.numRotors > 0:
            # Read the file names for the rotor logs
            line = readGeomFc.readMeaningfulLine(file)
            file_tokens = line.split()

            if len(file_tokens) != self.numRotors:
                print("Number of files doesn't match the number of rotors.")
                print("Exiting...")
                exit()

            # Read the symmetry factors for the internal rotors
            line = readGeomFc.readMeaningfulLine(file)
            sym_tokens = line.split()
            if (
                sym_tokens[0].upper() == "ROTORSYM"
                and len(sym_tokens) - 1 != self.numRotors
            ):
                print("Number of symmetry factors doesn't match the number of rotors.")
                print("Exiting...")
                exit()

            # Read the hindered rotor "frequency"
            # This will be removed from the frequencies list.
            line = readGeomFc.readMeaningfulLine(file)
            tokens = line.split()
            # print(len(tokens))

            if tokens[0].upper() == "ROTORFREQ" and len(tokens) - 1 != self.numRotors:
                print("Number of symmetry factors doesn't match the number of rotors.")
                print("Exiting...")
                exit()

            for i in range(self.numRotors):
                if self.verbose > 0:
                    # print(self.Freq)
                    print(
                        "Removing hindered rotor frequency %s cm^-1 from the list of vibrational frequencies."
                        % tokens[i + 1]
                    )
                rm_idx = self.Freq.index(float(tokens[i + 1]))
                self.hindFreq.append(self.Freq[rm_idx])  # Keep track of them
                del self.Freq[rm_idx]

            # Create rotor object for each hindered rotor
            self.rotors = []
            for i in range(self.numRotors):
                rotor_label = self.label + "_" + str(i + 1)
                self.rotors.append(
                    rotor.Rotor(
                        self.dir + file_tokens[i],
                        self.geom,
                        self.Mass,
                        self.bonds,
                        int(sym_tokens[i + 1]),
                        label=rotor_label,
                        v_ho=float(tokens[i + 1]),
                    )
                )

            # print(self.Freq) # TODO

            # exit(0)

    # read the bonds
    # line = readGeomFc.readMeaningfulLine(file)
    # tokens = line.split()
    # for bond in tokens:
    #     self.bonds.append(float(bond))

    ### end of read_input()

    def read_log(self, log_file):
        """
        Read all molecular data from single G16 log file.

        Arguments:

            log_file (string): Path to g16 log file. (Use absolute if possible).

        """

        with open(log_file, "r") as f:
            # Read geometry, mass, and bonds
            (self.geom, self.Mass, self.bonds) = readGeomFc.readGeom(f)
            self.calculateMomInertia()

        # TODO
        # Read frequencies
        freq_raw = readGeomFc.read_freq(log_file)
        if len(freq_raw) == 0 and self.Mass.size != 1:
            self.calc_complete = False
            return
        for freq in freq_raw:
            if freq < 0 and self.TS:
                self.imagFreq = freq
            elif freq < 0 and self.TS != True:
                self.calc_complete = False
                return
            else:
                self.Freq.append(freq)

        # TODO make this more flexible
        # Read energy
        self.Etype = "ub3lyp"
        self.Energy = readGeomFc.readEnergy(log_file, self.Etype)
        self.e_file = log_file

        # Dipole
        # TODO add doc. for this in class header
        self.dip_magn, self.dip = readGeomFc.read_dipole(log_file)  # in Debye

        # Misc
        self.extSymm = 1
        self.nelec = 2

        # TODO Add capability to handle rotors.
        self.numRotors = 0

    ############################################################################
    ### Primary Functions ######################################################
    ############################################################################
    def calculateMomInertia(self):
        geom = self.geom.copy()
        Mass = self.Mass.copy()
        geom /= 1e10
        Mass /= N_avo * 1e3

        # change coordinates to have center of mass
        cm = np.array([0, 0, 0])

        for i in range(Mass.size):
            cm = cm + Mass[i] * geom[i, :]

        cm = cm / np.sum(Mass)
        geom -= cm

        # calculate moments of inertia
        I = np.zeros((3, 3))
        x = geom[:, 0]
        y = geom[:, 1]
        z = geom[:, 2]

        I[0, 0] = np.sum(Mass * (y * y + z * z))
        I[1, 1] = np.sum(Mass * (x * x + z * z))
        I[2, 2] = np.sum(Mass * (x * x + y * y))
        I[0, 1] = I[1, 0] = -np.sum(Mass * x * y)
        I[0, 2] = I[2, 0] = -np.sum(Mass * x * z)
        I[1, 2] = I[2, 1] = -np.sum(Mass * z * y)

        # rotate coordinate axes to be parallel to principal axes
        (l, v) = la.eig(I)
        self.Iext = l  # in kg * m^2

    ############################################################################
    def calculate_thermo(self, oFile, Temp, mode_type="trans", harmonic=False):
        """
        Return the thermodynamic quantities at a given temperature for the type
        of mode specified.

        Arguments:

            oFile: open file object where the program will write the data

            Temp: list of temperatures in Kelvin

            type: (optional) the type of target mode ('trans', 'vib', 'rot',
                'int rot')

            harmonic: boolean flag for whether or not to calculate the
                properties harmonically for the internal rotor (only use if
                mode_type='int rot')

        Returns:

            Q: list of partition functions (for mode_type='trans' this is
                q/volume)

            H: list of enthalpies in kcal/mol (no ZPE included)

            Cp: list of heat capacities in cal/(mol K)

            S: list of entropies in cal/(mol K)

        """

        # Pint
        if oFile != False:
            if mode_type == "trans":
                self.print_thermo_heading(oFile, "Translational Contributions")
            elif mode_type == "vib":
                self.print_thermo_heading(oFile, "Vibrational Contributions")
            elif mode_type == "rot":
                self.print_thermo_heading(oFile, "Rotational Contributions")
            elif mode_type == "int rot":
                self.print_thermo_heading(oFile, "Int. Rotational Contributions")

        # Calculate the contributions at each temperature in the list
        Q = []
        H = []
        Cp = []
        S = []
        for T in Temp:
            # Ideal Gas Translational Modes
            # See Tuckerman, Stat. Mech.: Theory and Simulation section 5.6
            if mode_type == "trans":
                mass = np.sum(self.Mass) / 1e3 / N_avo  # Mass in kg / molecule
                Q_T = ((2 * np.pi * mass) / h ** 2) ** 1.5 / 101325 * (kb * T) ** 2.5
                H_T = 5.0 / 2.0 * R_kcal * T
                Cp_T = 5.0 / 2.0 * R_cal
                # See Sackur-Tetrode Equation
                S_T = R_cal * (np.log(Q_T) + 5.0 / 2.0)

            # Rigid Rotor Rotational Modes
            elif mode_type == "rot":
                sigma = self.extSymm
                Iext = self.Iext.copy()
                Iext[Iext == 0.0] = 1  # Prevents underflow if rotational symm
                # shows up TODO
                Q_T = np.power(np.pi * Iext[0] * Iext[1] * Iext[2], 0.5) / sigma
                Q_T *= np.power(8.0 * np.pi ** 2 * kb * T / h ** 2, 3.0 / 2.0)
                # print(Q_T)
                H_T = 3.0 / 2.0 * R_kcal * T
                Cp_T = 3.0 / 2.0 * R_cal
                # print(Iext, Q_T, sigma, Q_T/sigma) TODO
                S_T = H_T * 1e3 / T + R_cal * np.log(Q_T)

            # Harmonic Oscillator Vibrational Mode
            elif mode_type == "vib":
                freqs = np.array(self.Freq) * self.scale
                Q_T = 1.0
                H_T = 0
                Cp_T = 0
                S_T = 0
                for i in range(freqs.size):
                    ei = h * freqs[i] * c_in_cm  # hv for this mode in J
                    Q_T *= 1.0 / (1.0 - np.exp(-ei / (kb * T)))
                    H_T += ei / (np.exp(ei / (kb * T)) - 1.0) * (N_avo * j_to_cal / 1e3)
                    Cp_T += (
                        R_cal
                        * (ei / (kb * T)) ** 2
                        * np.exp(ei / (kb * T))
                        / (1.0 - np.exp(ei / (kb * T))) ** 2
                    )
                    S_T += R_cal * (
                        (ei / (kb * T)) / (np.exp(ei / (kb * T)) - 1.0)
                        - np.log(1.0 - np.exp(-ei / (kb * T)))
                    )

            # Internal Rotor Torsional Modes
            elif mode_type == "int rot":
                Q_T = 1
                H_T = 0
                Cp_T = 0
                S_T = 0

                for rotor in self.rotors:
                    q_ir, h_ir, cp_ir, s_ir = rotor.calculate_thermo(
                        T, harmonic=harmonic
                    )
                    Q_T *= q_ir
                    H_T += h_ir
                    Cp_T += cp_ir
                    S_T += s_ir

            Q.append(Q_T)
            H.append(H_T)
            Cp.append(Cp_T)
            S.append(S_T)
            # print("Cp %f" %Cp_T)
        # End of loop over temperatures
        # print(S)

        # Print thermo data
        if oFile != False:
            self.print_thermo_contributions(oFile, Temp, Q, H, Cp, S)
            oFile.write("\n")

        return Q, H, Cp, S

    ############################################################################
    def calculate_all_thermo(self, Temp, oFile):
        """
        Calculate and write the thermodynamic properties contributions from
        translation, vibration, rotation, and, if applicable, internal rotation.

        Arguments:

            Temp: list of temperatures

            oFile: an open file object where the thermo properties will be
                written

        Returns:

            Q: total partition function

            H: total enthalpy in kcal/mol (no ZPE included)

            Cp: total heat capacity in cal/(mol K)

            S: total entropy in cal/(mol K)
        """

        # Translation
        q_tr, h_tr, cp_tr, s_tr = self.calculate_thermo(oFile, Temp, mode_type="trans")

        # Rotation
        q_rot, h_rot, cp_rot, s_rot = self.calculate_thermo(
            oFile, Temp, mode_type="rot"
        )

        # Vibration
        q_vib, h_vib, cp_vib, s_vib = self.calculate_thermo(
            oFile, Temp, mode_type="vib"
        )

        # Internal Rotation
        q_ir, h_ir, cp_ir, s_ir = self.calculate_thermo(
            oFile, Temp, mode_type="int rot"
        )

        # Combine all of the data
        Q = [q_tr[i] * q_vib[i] * q_rot[i] * q_ir[i] for i in range(len(Temp))]
        H = [h_tr[i] + h_vib[i] + h_rot[i] + h_ir[i] for i in range(len(Temp))]
        Cp = [cp_tr[i] + cp_vib[i] + cp_rot[i] + cp_ir[i] for i in range(len(Temp))]
        S = [s_tr[i] + s_vib[i] + s_rot[i] + s_ir[i] for i in range(len(Temp))]

        if oFile != False:
            self.print_thermo_heading(oFile, "Total Thermo. Contributions")
            self.print_thermo_contributions(oFile, Temp, Q, H, Cp, S)

        return Q, H, Cp, S

    ############################################################################

    def print_thermo_heading(self, oFile, heading):
        """
        Write the header the thermo contributions of a particular set of modes
        (e.g. translationa, rotational, etc).

        Arguments:

            oFile: open output file object where program will write the header

            heading (string): string describing the thermo contributions to be
                written
        """

        l = len(heading) + 4
        symb = "-"
        header = footer = symb * l
        mainline = "\n%s %s %s\n" % (symb, heading, symb)
        oFile.write(header + mainline + footer + "\n\n")
        return

    ############################################################################

    def print_thermo_contributions(self, oFile, Temp, Q, H, Cp, S):
        """
        Write the thermo contributions from a set of modes to the cantherm
        output file.

        Arguments:

            oFile: open output file object

            Temp: list of temperatures for calculated contributions

            ent: entropy contributions

            cp: constant pressure heat capacity contributions

            dH: enthalpy contributions
        """

        horizontal_line = "-" * 12 + " " * 3
        horizontal_line *= 5
        horizontal_line += "\n"

        oFile.write(horizontal_line)
        temp_string = "Temp"
        temp_units_string = "K"
        q_string = "Q"
        q_units_string = "Unit-less"
        s_string = "S"
        s_units_string = "cal/(mol K)"
        cp_string = "Cp"
        cp_units_string = "cal/(mol K)"
        h_string = "H"
        h_units_string = "kcal/mol"

        oFile.write(
            "{:^12}   {:^12}   {:^12}   {:^12}   {:^12}\n".format(
                temp_string, q_string, h_string, cp_string, s_string
            )
        )
        oFile.write(
            "{:^12}   {:^12}   {:^12}   {:^12}    {:^12}\n".format(
                temp_units_string,
                q_units_string,
                h_units_string,
                cp_units_string,
                s_units_string,
            )
        )
        oFile.write(horizontal_line)

        for i in range(len(Temp)):
            oFile.write("{:^12.2f}".format(Temp[i]) + "   ")
            oFile.write("{:^12.6e}".format(Q[i]) + "   ")
            oFile.write("{:^12.6f}".format(H[i]) + "   ")
            oFile.write("{:^12.6f}".format(Cp[i]) + "   ")
            oFile.write("{:^12.6f}\n".format(S[i]))

        oFile.write(horizontal_line + "\n")
        return

    ############################################################################
    def print_properties(self, out_file):
        """
        Writes the data for the molecule after fully initializing molecule obj.
        """

        geom = self.geom
        Mass = self.Mass

        # Write coordinates with masses
        geom_heading = "Cartesian Geometry:\n"
        geom_heading += "-" * (len(geom_heading) - 1) + "\n\n"
        out_file.write(geom_heading)

        out_file.write(
            "%10s" % "Mass(amu)"
            + "%10s" % "X(ang)"
            + "%10s" % "Y(ang)"
            + "%10s" % "Z(ang)"
            + "\n"
        )
        for i in range(Mass.size):
            out_file.write(
                "%10.3f" % float(Mass[i])
                + "%10.4f" % float(geom[i, 0])
                + "%10.4f" % float(geom[i, 1])
                + "%10.4f" % float(geom[i, 2])
                + "\n"
            )
        out_file.write("\n\n")

        # Write vibrational frequencies
        vib_heading = "Vibrational Frequencies (cm-1):\n"
        vib_heading += "-" * (len(vib_heading) - 1) + "\n\n"
        out_file.write(vib_heading)

        Freq = self.Freq
        if self.TS:
            out_file.write("TS Imaginary Frequency: " + str(self.imagFreq) + "\n")

        for i in range(int(len(Freq) / 3) + 1):
            for j in range(3):
                if 3 * i + j < len(Freq):
                    out_file.write("%10.3f" % Freq[3 * i + j])
            out_file.write("\n")
        out_file.write("\n\n")

        # Write frequencies of hindered rotors excluded from q_vib
        hr_vib_heading = "Approx. Frequencies of Hindered Rotors (cm-1):\n"
        hr_vib_heading += "-" * (len(hr_vib_heading) - 1) + "\n\n"
        out_file.write(hr_vib_heading)

        for i in range(int(len(self.hindFreq) / 3) + 1):
            for j in range(3):
                if 3 * i + j < len(self.hindFreq):
                    out_file.write("%10.3f" % self.hindFreq[3 * i + j])
            out_file.write("\n")
        out_file.write("\n\n")

        # Write misc. properties
        misc_heading = "Other Properties:\n"
        misc_heading += "-" * (len(misc_heading) - 1) + "\n\n"
        out_file.write(misc_heading)

        if self.e_file != None:
            out_file.write("Energy file: %s\n" % self.e_file)
        out_file.write(
            "Energy = %10.3e kcal/mol %s\n" % (self.Energy * ha_to_kcal, self.Etype)
        )
        out_file.write("External Symmetry = " + str(self.extSymm) + "\n")
        Iext = self.Iext.copy() * 1e23 * N_avo
        out_file.write(
            "Principal Moments of Inertia (amu * ang.^2)= "
            + "%10.3f %10.3f %10.3f \n" % (Iext[0], Iext[1], Iext[2])
        )

        out_file.write("Electronic Degeneracy = " + str(self.nelec) + "\n\n")

        # Write Rotor properties
        for rotor in self.rotors:
            rotor.print_properties(out_file)

    ############################################################################
    ### Optional Functions #####################################################
    ############################################################################
    # def calculate_Q(self, T):
    #     """

    #     For more details see "Molecular Driving Forces" by Dill Chapters 11
    #     and 19.
    #     """

    #     # Translational Contrib.
    #     # TODO Assuming Unimolecular for now so it's technically per volume
    #     q_tr = np.power((2 * np.pi * sum(self.Mass) * kb * T) / h ** 2, 3.0 / 2.0)

    #     # Vibrational Contrib.
    #     q_vib = 1
    #     for freq in self.Freq:
    #         q_vib *= 1 / (1 - np.exp(-h * freq * c_in_cm / (kb * T)))

    #     # External Rotational Contrib.
    #     # TODO Add symm. factor
    #     sigma = 1
    #     q_rot = (
    #         np.power(np.pi * self.Iext[0] * self.Iext[1] * self.Iext[2], 0.5) / sigma
    #     )
    #     q_rot *= np.power(8.0 * np.pi ** 2 * kb * T / h ** 2, 3.0 / 2.0)

    #     # Internal Rotational Contrib.
    #     q_ir = 1
    #     s_ir = 0
    #     h_ir = 0
    #     cp_ir = 0

    #     for rotor in self.rotors:

    #         q_ir_i, s_ir_i, h_ir_i, cp_ir_i = rotor.calculate_thermo(T)

    #         q_ir *= q_ir_i
    #         s_ir += s_ir_i
    #         h_ir += s_ir_i
    #         cp_ir += s_ir_i

    #     return q_tr * q_vib * q_rot * q_ir
