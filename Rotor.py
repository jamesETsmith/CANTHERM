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
import numpy as np
from functools import reduce
la = np.linalg

import readGeomFc

def er_rotation(v1, v2):
    '''
    Euler-Rodrigues rotation of vector 1 to align with vector 2.

    Arguments:
        v1: vector that will be rotated
        v2: vector that we will rotate to (i.e. we will make v1 || to v2)

    Returns:
        r: 3x3 rotation matrix
    '''

    # Vector we will rotate about
    k = np.cross(v1,v2)
    # Angle we need to rotate
    th = np.arccos(np.dot(v1,v2)/(la.norm(v1)*la.norm(v2)))

    # Euler/Rodrigues params
    # See https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula
    a = np.cos(th/2)
    b = k[0] * np.sin(th/2)
    c = k[1] * np.sin(th/2)
    d = k[2] * np.sin(th/2)

    r = np.zeros((3,3))
    r[0,0]=a**2 + b**2 - c**2 - d**2; r[0,1]=2*(b*c-a*d); r[0,2]=2*(b*d+a*c);
    r[1,0]=2*(b*c+a*d); r[1,1]=a**2 + c**2 - b**2 - d**2; r[1,2]=2*(c*d-a*b);
    r[2,0]=2*(b*d-a*c); r[2,1]=2*(c*d+a*b); r[2,2]=a**2 + d**2 - b**2 - c**2;

    return r

class Rotor:
    '''
    Class for the hindered rotors.

       .. note:

           Make sure the atoms are in the same order as the original energy
    '''
    pivotAtom = 0
    atomsList = []
    pivot2 = 0
    moments = []
    level = 0
    r = mat('0.0 0.0 0.0')
    symm = 1.0
    energies = []
    pot_coeffs = np.array([])

    def __init__(self, log_file, geom, masses, bonds):
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

        self.geom = geom # list
        self.masses = masses  # 2D list index 1 is atom, index 2 is x/y/z
        self.natom = len(self.masses)
        self.bonds = bonds

        self.log_file =  log_file
        self.read_scan_data(self.log_file)
        self.fit_potential()
        # self.plot()
        self.calc_axis()
        # self.plot_rotational_pes()

    # def getAxes(self, geom, Mass):
    #     z = -(geom[self.pivot2 - 1, :] - geom[self.pivotAtom - 1, :])
    #     z = z / linalg.norm(z)
    #
    #     cm = matrix('0.0 0.0 0.0')
    #
    #     M = 0.0
    #     for i in self.atomsList:
    #         cm = cm + Mass[i - 1] * geom[i - 1, :]
    #         M = M + Mass[i - 1]
    #     cm = cm / M
    #
    #     xtemp = (cm - geom[self.pivotAtom - 1, :])
    #     xtemp = xtemp / linalg.norm(xtemp)
    #
    #     diff = xtemp - z
    #     different = False
    #
    #     for i in range(3):
    #         if not(-1e-10 < (xtemp[0, i] - z[0, i]) < 1e-10):
    #             different = True
    #             break
    #
    #     if (different):
    #         x = xtemp - (z * transpose(xtemp)) * z
    #         x = x / linalg.norm(x)
    #         y = matrix(cross(z, x))
    #     else:
    #         xtemp = z + mat(' 0.0 0.0 1.0')
    #         x = xtemp - (z * transpose(xtemp)) * z
    #         x = x / linalg.norm(x)
    #         y = matrix(cross(z, x))
    #
    #     self.dircos = matrix(zeros((3, 3), dtype=float))
    #     self.dircos[0, :] = x
    #     self.dircos[1, :] = y
    #     self.dircos[2, :] = z

    def calc_axis(self):
        '''
        Calculates the axis of hindered rotation and assign mocules to their
        respective top (i.e. which side of the bond they are on).
        '''

        geom = np.array(self.geom)
        masses = np.array(self.masses)
        natom = self.natom

        # Option 1: Along the axis of the dihedral bond
        # Shift pivot 2 to the origin
        p2 = np.array(geom[self.pivot2,:])
        geom[:] -= p2
        p1 = np.array(geom[self.pivot1,:]) # After translation

        rot_mat = er_rotation(p1,np.array([0,0,1])) # Rotate to z-axis
        self.geom = np.dot(geom,rot_mat.T) #np.einsum('ij,kj->ik', rot_mat, geom)
        # self.plot() # TODO

        # Assign the atoms to atomlist for each top
        self.alist1 = [] # closest to pivot 1
        self.alist2 = [] # closest to pivot 2

        sets = [ [i] for i in range(natom) ]
        # print("SETS0", sets) # TODO
        idx = self.bonds.index([min(self.pivot1,self.pivot2),
                                max(self.pivot1,self.pivot2)])
        del self.bonds[idx]

        for b in self.bonds:
            idx0 = None
            idx1 = None

            for i in range(len(sets)):
                if b[0] in sets[i]:
                    idx0 = i
                if b[1] in sets[i]:
                    idx1 = i

            # Already in the same set
            if idx1 == idx0:
                continue

            # Combine the sets and clean up
            else:
                sets[idx0] += sets[idx1]
                del sets[idx1]

        # print("SETS", sets) # TODO
        self.bonds.insert( idx, [min(self.pivot1,self.pivot2),
                                max(self.pivot1,self.pivot2)])

        # Collect the list of atoms in each top
        p1_idx = None
        p2_idx = None
        for i in range(len(sets)):
            if self.pivot1 in sets[i]:
                p1_idx = i
            if self.pivot2 in sets[i]:
                p2_idx = i

        self.alist1 = sets[p1_idx]
        self.alist2 = sets[p2_idx]
        self.alist1.remove(self.pivot1)
        self.alist2.remove(self.pivot2)

        print("ALIST1", self.alist1, self.pivot1, len(self.alist1))
        print("ALIST2", self.alist2, self.pivot2)
        return

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
                if re.search('D\s+\d', line):
                    print(line.split()) # TODO
                    splt = line.split()
                    self.dihedral = np.array([int(splt[1]), int(splt[2]),
                                              int(splt[3]), int(splt[4])]) - 1
                    self.pivot1 = self.dihedral[1]
                    self.pivot2 = self.dihedral[2]
                    self.scan_rate = float(splt[-1]) # In degrees

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
        return


    def plot(self):
        '''
        Plot the rotors and highlight them.

        '''

        from mpl_toolkits.mplot3d import axes3d

        geom = self.geom
        masses = self.masses

        # Figure
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Add Atoms
        coords = np.array(geom)
        # print(coords[:,0]) #TODO

        for i in range(len(masses)):
            if masses[i] == 12.0: # Carbon grey
                ax.scatter(coords[i,0], coords[i,1], coords[i,2], c='#C8C8C8',
                    s=300, depthshade=False)
            elif masses[i] == 15.99491: # Oxygen red
                ax.scatter(coords[i,0], coords[i,1], coords[i,2], c='#F00000',
                    s=300, depthshade=False)
            elif masses[i] == 1.00783: # Hydrogen white
                ax.scatter(coords[i,0], coords[i,1], coords[i,2], c='#FFFFFF',
                    s=300, depthshade=False)
            elif masses[i] == 14.0031: # Nitrogen light blue
                ax.scatter(coords[i,0], coords[i,1], coords[i,2], c='#8F8FFF',
                    s=300, depthshade=False)
            elif masses[i] == 34.96885: # Chlorine green
                ax.scatter(coords[i,0], coords[i,1], coords[i,2], c='#00FF00',
                    s=300, depthshade=False)
            elif masses[i] == 31.97207: # Sulphur yellow
                ax.scatter(coords[i,0], coords[i,1], coords[i,2], c='#FFC832',
                    s=300, depthshade=False)
            elif masses[i] == 18.99840: # Fluorine seafoam green
                ax.scatter(coords[i,0], coords[i,1], coords[i,2], c='#0fff97',
                    s=300, depthshade=False)

        ax.scatter(coords[self.pivot1,0], coords[self.pivot1,1],
                   coords[self.pivot1,2], c='#6f6f6f', s=300, depthshade=False,
                   marker = '^')

        ax.scatter(coords[self.pivot2,0], coords[self.pivot2,1],
                   coords[self.pivot2,2], c='#6f6f6f', s=300, depthshade=False,
                   marker = '^')

        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')

        plt.show()

        return


    def write_potential_to_file(self, energy_offset, energies, file_name):
        pot_file_header = "# 1-D Scan of Total Energy\n#\n# %8s\t%16s\n"%('Step Number','Energy (Ha)')
        with open(file_name, 'w') as f:
            f.write(pot_file_header)
            i = 1;
            for energy in energies:
                f.write("%8i\t\t\t%16.8f\n"%(i,energy+energy_offset))
                i += 1

        return

    def plot_rotational_pes(self, fitted=True, n_fit=100):
        '''
        Plot the 1-D Rotor PES.
        '''

        if fitted:
            th = np.linspace(0, 2*np.pi, num=n_fit)
            v_fit = self.calculate_potential(th=th)

        plt.figure()
        plt.plot(np.linspace(0,2*np.pi, num=len(self.energies)),
                 np.array(self.energies), # - np.min(self.energies),# + self.energy_offset,
                 label="Calculated")
        plt.plot(th, v_fit, label="Fit")
        plt.legend()
        plt.show()


    def fit_potential(self, a=None, e=None, n_term=23):
        '''
        A is the matrix of fourier terms. E are the energies. N_term is the
        number of Fourier terms to include in the expansion.
        '''

        if e == None:
            e = np.array(self.energies) #- np.min(self.energies) # self.energy_offset
            e = e[:-1] # Last element is a duplicate.

        # e = np.sin(np.linspace(0,2*np.pi,num=24))/100  + np.sin(3 * np.linspace(0,2*np.pi,num=24))/100
        # self.energies = e

        if a == None:
            a = np.zeros( (e.shape[0], n_term) )
            th = np.linspace(0,2*np.pi,e.shape[0])

            for r in range(a.shape[0]):
                for c in range(a.shape[1]):
                    # First coefficient (no sinusoidal function)
                    if c == 0:
                        a[r,c] = 1
                    # a coeffs
                    elif c%2 == 0 and c > 0:
                    #    a[r,c] = 1 - np.cos(k[c]/2 * th[r])
                        a[r,c] = np.cos(int((c + 1)/2) * th[r])
                    # b coeffs
                    else:
                        a[r,c] = np.sin(int((c + 1)/2) * th[r])

        self.pot_coeffs = reduce(np.dot, (la.inv(np.dot(a.T, a)), a.T, e))
        # print(self.pot_coeffs) # TODO
        return

    def calculate_potential(self, th=None, n_fit=1000 ):
        '''
        Calculate the fitted potential using the pot_coeffs array.
        '''

        if th is None:
            th = np.linspace(0, 2*np.pi, num=n_fit)

        v_fit = np.zeros(th.shape[0])

        if self.pot_coeffs.shape[0] == 0: self.fit_potential()

        for i in range(th.shape[0]):
          for j in range(self.pot_coeffs.shape[0]):
              if j == 0:
                  v_fit[i] += self.pot_coeffs[j] * 1.0
                  continue
              elif j%2 == 0 and j > 0:
                  v_fit[i] += self.pot_coeffs[j] * \
                      np.cos( int((j + 1)/2) * th[i])
              elif j%2 == 1:
                  v_fit[i] += self.pot_coeffs[j] * \
                      np.sin( int((j + 1)/2) * th[i])

        return v_fit
