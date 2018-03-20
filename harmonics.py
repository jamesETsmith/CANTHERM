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

from numpy import *
import pdb
import readGeomFc
import math


class Harmonics:
    numFit = 0
    Kcos = []
    Ksin = []
    A = 0.0

    def __init__(self, numFit, Kcos, Ksin):
        self.numFit = numFit
        self.Kcos = Kcos
        self.Ksin = Ksin

    def getPotential(self, angle):
        # assuming the angle is in degrees
        rad = angle * math.pi / 180
        pot = self.A
        for i in range(self.numFit):
            pot = pot + \
                self.Ksin[i] * math.sin((i + 1) * rad) + \
                self.Kcos[i] * math.cos((i + 1) * rad)
        return pot

    def fitPotential(self, file):
        read = open(file, 'r')
        lines = read.readlines()
        angles = []
        potentials = []
        pot0 = 0.0
        # first three lines are comments and last line is repeat
        nfit = len(lines) - 4
        potgiven = []
        for i in range(3, len(lines) - 1):
            tokens = lines[i].split()

            potentials.append(float(tokens[1]))
            if i == 3:
                pot0 = potentials[0]
            potentials[i - 3] = (potentials[i - 3] - pot0) * ha_to_kcal
            potgiven.append([float(i - 3) * 360 / nfit, potentials[i - 3]])

        # now fit the potentials
        Y = transpose(matrix(potentials[:nfit]))
        X = matrix(zeros((nfit, 11), dtype=float))
        for i in range(nfit):
            angle = float(i) * 2 * math.pi / nfit
            X[i, 0] = 1.0
            for j in range(5):
                X[i, j + 1] = math.cos((j + 1) * angle)
                X[i, 6 + j] = math.sin((j + 1) * angle)

        XtX = transpose(X) * X
        XtY = transpose(X) * Y
        b = linalg.inv(XtX) * XtY

        for i in range(5):
            self.Kcos.append(0.0)
            self.Kcos[i] = float(b[i + 1])
            self.Ksin.append(0.0)
            self.Ksin[i] = float(b[i + 6])
        self.A = float(b[0])
        self.numFit = 5

        pot = []
        for i in range(3 * nfit + 1):
            angle = i * 360 / 3 / nfit
            pot.append([angle, self.getPotential(angle)])

        return
