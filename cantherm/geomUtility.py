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
from harmonics import *
from scipy import *
import random
import math


def rotateDihedrals(geom, rotors, dihedrals, Mass):
    nGeom = geom.copy()
    # print geom
    k = 0
    for rotor in rotors[1:]:
        # move center to pivot atom
        pivotcoord = nGeom[rotor.pivotAtom - 1, :].copy()
        for i in range(Mass.size):
            nGeom[i, :] = nGeom[i, :] - pivotcoord

    # make the axis of rotor the z axis
        # first generate the required direction cosines of the new coordinates
        rotor.getAxes(nGeom, Mass)
        for i in range(Mass.size):
            nGeom[i, :] = transpose(rotor.dircos * transpose(nGeom[i, :]))

    # for this rotor rotate the atoms in x,y plane by dihedrals[i] degrees
        dihed = dihedrals[k]
        for atoms in rotor.atomsList[1:]:
            x = nGeom[atoms - 1, 0]
            y = nGeom[atoms - 1, 1]
            r = sqrt(x**2 + y**2)
            theta = math.acos(x / r)
            if y < 0.0:
                theta = -theta
            theta = theta + dihed * math.pi / 180.0
            x = r * math.cos(theta)
            y = r * math.sin(theta)
            nGeom[atoms - 1, 0] = x
            nGeom[atoms - 1, 1] = y
        k = k + 1
    return nGeom

#***********************************************************************************


def findCM(geom, Mass):
    cm = matrix('0.0 0.0 0.0')

    for i in range(Mass.size):
        cm = cm + Mass[i] * geom[i, :]

    cm = cm / sum(Mass)

    return cm


#***********************************************************************************

def calculateD32(geom, Mass, rotors):
    numRotors = len(rotors)

    redMom = matrix(zeros((numRotors - 1, 1), dtype=float))
    k = 0
    for rotor in rotors[1:]:
        RotorList = rotor.atomsList
        nonRotorList = []
        for j in range(len(Mass)):
            if (RotorList.count(j + 1) == 0):
                nonRotorList.append(j + 1)

        # calculate the cm for rotor and nonRotorList
        cm1 = mat('0.0 0.0 0.0')
        mass1 = 0.0
        cm2 = mat('0.0 0.0 0.0')
        mass2 = 0.0
        for j in RotorList:
            cm1 = cm1 + Mass[j - 1] * geom[j - 1, :]
            mass1 += Mass[j - 1]
        cm1 = cm1 / mass1

        for j in nonRotorList:
            cm2 = cm2 + Mass[j - 1] * geom[j - 1, :]
            mass2 += Mass[j - 1]
        cm2 = cm2 / mass2

        axofrot = (cm1 - cm2) / linalg.norm(cm1 - cm2)  # axis of rotation

        I1 = 0.0  # moment of inertia of first fragment
        I2 = 0.0  # moment of inertia of second fragment

        for j in RotorList:
            r1 = (geom[j - 1, :] - cm1) - \
                ((geom[j - 1, :] - cm1) * transpose(axofrot)) * axofrot
            I1 += Mass[j - 1] * linalg.norm(r1)**2

        for j in nonRotorList:
            r1 = (geom[j - 1, :] - cm2) - \
                ((geom[j - 1, :] - cm2) * transpose(axofrot)) * axofrot
            I2 += Mass[j - 1] * linalg.norm(r1)**2

        redMom[k] = 1.0 / (1.0 / I1 + 1.0 / I2)
        k += 1
    return redMom

#***********************************************************************************


def calculateD(geom, Mass, rotors):
    numRotors = len(rotors)

# change coordinates to have cm
    cm = matrix('0.0 0.0 0.0')

    for i in range(Mass.size):
        cm = cm + Mass[i] * geom[i, :]

    cm = cm / sum(Mass)

    for i in range(Mass.size):
        geom[i, :] = geom[i, :] - cm


# calculate moments of inertia
    I = matrix(zeros((3, 3), dtype=double))
    x = array(geom[:, 0])
    y = array(geom[:, 1])
    z = array(geom[:, 2])
    I[0, 0] = sum(array(Mass) * (y * y + z * z))
    I[1, 1] = sum(array(Mass) * (x * x + z * z))
    I[2, 2] = sum(array(Mass) * (x * x + y * y))
    I[0, 1] = I[1, 0] = -sum(array(Mass) * x * y)
    I[0, 2] = I[2, 0] = -sum(array(Mass) * x * z)
    I[1, 2] = I[2, 1] = -sum(array(Mass) * z * y)

# rotate coordinate axes to be parallel to principal axes
    (l, v) = linalg.eig(I)
    prinAxes = transpose(v)

    for i in range(Mass.size):
        geom[i, :] = transpose(prinAxes * transpose(geom[i, :]))

# again calculate the moments of inertia and confirm that they are diagonal
    I = matrix(zeros((3, 3), dtype=double))
    x = array(geom[:, 0])
    y = array(geom[:, 1])
    z = array(geom[:, 2])
    I[0, 0] = sum(array(Mass) * (y * y + z * z))
    I[1, 1] = sum(array(Mass) * (x * x + z * z))
    I[2, 2] = sum(array(Mass) * (x * x + y * y))
    I[0, 1] = I[1, 0] = -sum(array(Mass) * x * y)
    I[0, 2] = I[2, 0] = -sum(array(Mass) * x * z)
    I[1, 2] = I[2, 1] = -sum(array(Mass) * z * y)
#    print I

    K = matrix(zeros((6 + numRotors - 1, 6 + numRotors - 1), dtype=double))

    M = eye(3, 3)
    M = matrix(M * sum(Mass))
    U1 = matrix('0.0; 0.0; 0.0')

    for i in range(3):
        U1[i] = sum(array(Mass) * array(geom[:, i]))

    U = matrix(zeros((3, 3), dtype=double))
    U[0, 1] = -1 * U1[2]
    U[1, 0] = U1[2]
    U[0, 2] = U1[1]
    U[2, 0] = -1 * U1[1]
    U[1, 2] = -1 * U1[0]
    U[2, 1] = U1[0]

    K[0:3, 0:3] = M
    K[3:6, 0:3] = U
    K[0:3, 3:6] = transpose(U)
    K[3:6, 3:6] = I

    irotor = 0
    for rotor in rotors[1:]:
        rotor.getAxes(geom, Mass)
        rotor.getMoments(geom, Mass)

        A = rotor.moments[0]
        B = rotor.moments[1]
        C = rotor.moments[2]
        Ux = rotor.moments[3]
        dircos = rotor.dircos
        r = transpose(rotor.r)
        level = rotor.level

    #    print dircos[2,:]
        K[6 + irotor, 0:3] = Ux * dircos[1, :]
        K[0:3, 6 + irotor] = transpose(Ux * dircos[1, :])

        beta_i1 = mat('0.0;0.0;0.0')
        for i in range(3):
            beta_i1[i] = dircos[2, i] * A - dircos[0, i] * B - dircos[1, i] * C

        beta_i1[0] = beta_i1[0] + \
            (dircos[1, 2] * r[1] - dircos[1, 1] * r[2]) * Ux
        beta_i1[1] = beta_i1[1] + \
            (dircos[1, 0] * r[2] - dircos[1, 2] * r[0]) * Ux
        beta_i1[2] = beta_i1[2] + \
            (dircos[1, 1] * r[0] - dircos[1, 0] * r[1]) * Ux

        K[6 + irotor, 3:6] = transpose(beta_i1)
        K[3:6, 6 + irotor] = beta_i1

        K[6 + irotor, 6 + irotor] = A

        if (level >= 3):
            pivot2 = rotor.pivot2
            numAncestors = level - 2

            presentParent = rotor
            # print numAncestors
            for i in range(numAncestors):
                parentNum = presentParent.parent

                presentParent = rotors[parentNum]
                dircos_ipar = dircos * transpose(presentParent.dircos)

                ripar = transpose((geom[rotor.pivotAtom - 1, :] -
                                   geom[presentParent.pivotAtom - 1, :]) * presentParent.dircos)

                beta_z = dircos_ipar[2, 2] * A - \
                    dircos_ipar[0, 2] * B - dircos_ipar[1, 2] * C
                beta_z = beta_z + (dircos_ipar[1, 1] *
                                   ripar[0] - dircos[1, 0] * ripar[1]) * Ux
                K[6 + irotor, 6 + parentNum - 1] = K[6 +
                                                     parentNum - 1, 6 + irotor] = beta_z

        irotor = irotor + 1

    S = K[3:, 3:] - K[3:, 0:3] * linalg.inv(K[0:3, 0:3]) * K[0:3, 3:]
    D = S[3:, 3:] - S[3:, 0:3] * linalg.inv(S[0:3, 0:3]) * S[0:3, 3:]
    return D


#***********************************************************************************


def calculateI43(geom, Mass, rotors):
    numRotors = len(rotors)

# change coordinates to have cm
    cm = matrix('0.0 0.0 0.0')

    for i in range(Mass.size):
        cm = cm + Mass[i] * geom[i, :]

    cm = cm / sum(Mass)

    for i in range(Mass.size):
        geom[i, :] = geom[i, :] - cm


# calculate moments of inertia
    I = matrix(zeros((3, 3), dtype=double))
    x = array(geom[:, 0])
    y = array(geom[:, 1])
    z = array(geom[:, 2])
    I[0, 0] = sum(array(Mass) * (y * y + z * z))
    I[1, 1] = sum(array(Mass) * (x * x + z * z))
    I[2, 2] = sum(array(Mass) * (x * x + y * y))
    I[0, 1] = I[1, 0] = -sum(array(Mass) * x * y)
    I[0, 2] = I[2, 0] = -sum(array(Mass) * x * z)
    I[1, 2] = I[2, 1] = -sum(array(Mass) * z * y)

# rotate coordinate axes to be parallel to principal axes
    (l, v) = linalg.eig(I)
    prinAxes = transpose(v)

    for i in range(Mass.size):
        geom[i, :] = transpose(prinAxes * transpose(geom[i, :]))

# again calculate the moments of inertia and confirm that they are diagonal
    I = matrix(zeros((3, 3), dtype=double))
    x = array(geom[:, 0])
    y = array(geom[:, 1])
    z = array(geom[:, 2])
    I[0, 0] = sum(array(Mass) * (y * y + z * z))
    I[1, 1] = sum(array(Mass) * (x * x + z * z))
    I[2, 2] = sum(array(Mass) * (x * x + y * y))
    I[0, 1] = I[1, 0] = -sum(array(Mass) * x * y)
    I[0, 2] = I[2, 0] = -sum(array(Mass) * x * z)
    I[1, 2] = I[2, 1] = -sum(array(Mass) * z * y)
#    print I

    K = matrix(zeros((6 + numRotors - 1, 6 + numRotors - 1), dtype=double))

    M = eye(3, 3)
    M = matrix(M * sum(Mass))
    U1 = matrix('0.0; 0.0; 0.0')

    for i in range(3):
        U1[i] = sum(array(Mass) * array(geom[:, i]))

    U = matrix(zeros((3, 3), dtype=double))
    U[0, 1] = -1 * U1[2]
    U[1, 0] = U1[2]
    U[0, 2] = U1[1]
    U[2, 0] = -1 * U1[1]
    U[1, 2] = -1 * U1[0]
    U[2, 1] = U1[0]

    K[0:3, 0:3] = M
    K[3:6, 0:3] = U
    K[0:3, 3:6] = transpose(U)
    K[3:6, 3:6] = I

    irotor = 0
    beta_i = []
    redMom = matrix(zeros((numRotors - 1, numRotors - 1), dtype=double))
    Amm = matrix(zeros((numRotors - 1, numRotors - 1), dtype=double))
    Im0 = mat(zeros((numRotors - 1, 1), dtype=double))

    for rotor in rotors[1:]:
        rotor.getAxes(geom, Mass)
        rotor.getMoments(geom, Mass)

        A = rotor.moments[0]
        B = rotor.moments[1]
        C = rotor.moments[2]
        Ux = rotor.moments[3]
        dircos = rotor.dircos
        r = transpose(rotor.r)

        beta_i1 = mat('0.0;0.0;0.0')
        for i in range(3):
            beta_i1[i] = dircos[2, i] * A - dircos[0, i] * B - dircos[1, i] * C

        beta_i1[0] = beta_i1[0] + \
            (dircos[1, 2] * r[1] - dircos[1, 1] * r[2]) * Ux
        beta_i1[1] = beta_i1[1] + \
            (dircos[1, 0] * r[2] - dircos[1, 2] * r[0]) * Ux
        beta_i1[2] = beta_i1[2] + \
            (dircos[1, 1] * r[0] - dircos[1, 0] * r[1]) * Ux

        beta_i.append(beta_i1)

        Im0[irotor] = A
        for i in range(3):
            Im0[irotor] = Im0[irotor] - \
                (dircos[1, i] * Ux)**2 / sum(Mass) - beta_i1[i]**2 / I[i, i]
        irotor = irotor + 1

    for k in range(len(rotors) - 1):
        for l in range(len(rotors) - 1):
            rotk = rotors[k + 1]
            rotl = rotors[l + 1]
            for i in range(3):
                Amm[k, l] = Amm[k, l] + rotk.dircos[1, i] * rotl.dircos[1, i] * rotk.moments[3] * rotl.moments[3] / sum(Mass) \
                    + beta_i[k][i] * beta_i[l][i] / I[i, i]

    for i in range(len(rotors) - 1):
        redMom[i, i] = Im0[i]
    return redMom
