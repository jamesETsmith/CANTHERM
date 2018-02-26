#!/usr/bin/env python
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

from numpy import *
import pdb
import re #, os, sys
import matplotlib.pyplot as plt



class Rotor:

    pivotAtom = 0
    atomsList = []
    pivot2 = 0
    moments = []
    level = 0
    r = mat('0.0 0.0 0.0')
    symm = 1.0
    energies = []

    def __init__(self, log_file):
        # self.pivotAtom = atomsList[0]
        # self.atomsList = atomsList
        # self.pivot2 = pivot2
        # self.level = level
        # self.parent = 0
        # self.symm = symm
        # self.nonRotorList = []
        # for j in range(len(Mass)):
        #     if (self.atomsList.count(j + 1) == 0):
        #         self.nonRotorList.append(j + 1)
        self.log_file =  log_file
        self.read_scan_data(self.log_file)
        self.plot_rotational_pes()

    def getAxes(self, geom, Mass):
        z = -(geom[self.pivot2 - 1, :] - geom[self.pivotAtom - 1, :])
        z = z / linalg.norm(z)

        cm = matrix('0.0 0.0 0.0')

        M = 0.0
        for i in self.atomsList:
            cm = cm + Mass[i - 1] * geom[i - 1, :]
            M = M + Mass[i - 1]
        cm = cm / M

        xtemp = (cm - geom[self.pivotAtom - 1, :])
        xtemp = xtemp / linalg.norm(xtemp)

        diff = xtemp - z
        different = False

        for i in range(3):
            if not(-1e-10 < (xtemp[0, i] - z[0, i]) < 1e-10):
                different = True
                break

        if (different):
            x = xtemp - (z * transpose(xtemp)) * z
            x = x / linalg.norm(x)
            y = matrix(cross(z, x))
        else:
            xtemp = z + mat(' 0.0 0.0 1.0')
            x = xtemp - (z * transpose(xtemp)) * z
            x = x / linalg.norm(x)
            y = matrix(cross(z, x))

        self.dircos = matrix(zeros((3, 3), dtype=float))
        self.dircos[0, :] = x
        self.dircos[1, :] = y
        self.dircos[2, :] = z

    def getMoments(self, geom, Mass):
        geomTemp = geom.copy()

        r = geom[self.pivotAtom - 1, :]

        self.r = r.copy()

        # translate so that pivot atom is at the origin
        for i in range(Mass.size):
            geomTemp[i, :] = geom[i, :] - geom[self.pivotAtom - 1, :]

        # now rotate so that the axes are parallel to rotor axes
        for i in range(Mass.size):
            geomTemp[i, :] = transpose(self.dircos * transpose(geomTemp[i, :]))

        A = 0.0  # moment of inertia about the z axis
        B = 0.0  # xz cross product of ineria
        C = 0.0  # yz cross product of inertia
        Ux = 0.0  # first order x moment
        x = geomTemp[:, 0]
        y = geomTemp[:, 1]
        z = geomTemp[:, 2]

        Uy = 0.0

        for k in self.atomsList:
            i = k - 1
            A = A + Mass[i] * (x[i]**2 + y[i]**2)
            B = B + Mass[i] * x[i] * z[i]
            C = C + Mass[i] * y[i] * z[i]
            Ux = Ux + Mass[i] * x[i]
            Uy = Uy + Mass[i] * y[i]

        self.moments = [A, B, C, Ux]

    def read_scan_data(self, file_name):
        energy_offset = 0
        energies = []

        # Step through file one line at a time
        with open(file_name, 'r') as f:
            lines = f.readlines()
            parsing_rotor = False #Boolean for parsing PES energy data


            for line in lines:
                # Get pivot atoms
                if re.match('D\s+', line):
                    print(line.split())

                # Start of scan data
                if re.match(' Summary of Optimized Potential Surface Scan', line):
                    energy_offset = float(line.split()[7])
                    parsing_rotor = True

                # Energies
                if parsing_rotor:
                    if re.match('\s+Eigenvalues --', line):
                        for e in line.split()[2:]:
                            energies.append(float(e))

                # End of scan data
                if parsing_rotor and re.match(' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad', line):
                    parsing_rotor = False
                    break

        # TODO Clean this up and just use self.* throughout function
        self.energies = energies
        self.energy_offset = energy_offset

    def write_potential_to_file(self, energy_offset, energies, file_name):
        pot_file_header = "# 1-D Scan of Total Energy\n#\n# %8s\t%16s\n"%('Step Number','Energy (Ha)')
        with open(file_name, 'w') as f:
            f.write(pot_file_header)
            i = 1;
            for energy in energies:
                f.write("%8i\t\t\t%16.8f\n"%(i,energy+energy_offset))
                i += 1

        return

    def plot_rotational_pes(self):
        '''
        Plot the 1-D Rotor PES.
        '''

        plt.figure()
        plt.plot(self.energies)
        plt.show()
