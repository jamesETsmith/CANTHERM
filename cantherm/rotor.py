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

### General Modules ###
import re
import numpy as np
from functools import reduce
import matplotlib.pyplot as plt
la = np.linalg

### CanTherm Modules ###
from cantherm.constants import *

class Rotor:
    '''
    Attributes:

        geom (np.array): 2D numpy array of floats that holds the coords in
            angstroms.

        mass (np.array): 1D numpy array of floats for the masses, index matches
            atom number.

        natom (int): Number of atoms in molecule.

        bonds ([[int]]): 2D list of the bonds (atoms are the atom number - 1).

        log_file (string): Path to scan file for this rotor. (Not file object).

        sym (float): Symmetry of hindered rotor.

        v_ho (float): Harmonic oscillator frequency associated with rotor.

        dihedral ([int]): List of atom indices used in dihedral scan.

        pivot1 (int): Atom index for the first pivot of internal rotor.

        pivot2 (int): Atom index for the second pivot of internal rotor.

        scan_rate (float): Number of degrees scanned per step.

        energy_offset (float): Energy offset in J for the rotational PES
            energies.

        energies (np.array): 1D numpy array of floats. Full energies at each
            step of the rotational PES in J.

    Saved Results:

        v0_coeff (float): Fourier coefficient for non trig. parameter in J.

        pot_coeffs (np.array): 2D numpy array of floats. Fourier coefficients
            for the cos(kx) and sin(kx) terms. pot_coeff[0] are the cosine terms
            and pot_coeff[1] are the sine terms. All coefficients are in units
            of J.

        i_red (float): Reduced moment of inertia for the rotor in kg*m^2.

        e_levels (np.array): 1D numpy array of floats. Energies of internal
            rotor (eigenvalues of the interal rotor hamiltonian).

    .. note:

        Make sure the atoms are in the same order as the original energy
    '''

    def __init__(self, log_file, geom, masses, bonds, sym, label='', v_ho=0):
        # Member variables
        self.geom = geom # list
        self.masses = masses #np.array( [masses[i,0] for i in range(len(masses))] )
        self.natom = self.masses.size
        self.bonds = bonds
        self.log_file =  log_file
        self.sym = sym
        self.e_levels = np.array([])
        self.v_ho = v_ho
        self.label = label

        # Functions
        self.read_scan_data(self.log_file)
        self.fit_potential()
        # self.plot()
        self.load_atom_lists()
        self.plot_rotational_pes()
        self.calc_reduced_moment()
        self.calc_e_levels()

    ######################################################################
    ### Mandatory Functions ##############################################
    ######################################################################
    def load_atom_lists(self):
        '''
        Use bonds to assign which atoms are a part of the rotor.
        '''

        # geom = np.array(self.geom)
        natom = self.natom

        # Option 1: Along the axis of the dihedral bond
        # Shift pivot 2 to the origin
        # p2 = np.array(geom[self.pivot2,:])
        # geom[:] -= p2
        # p1 = np.array(geom[self.pivot1,:]) # After translation
        #
        # rot_mat = er_rotation(p1,np.array([0,0,1])) # Rotate to z-axis
        # self.geom = np.dot(geom,rot_mat.T) #np.einsum('ij,kj->ik', rot_mat, geom)
        # self.plot() # TODO

        # Assign the atoms to atomlist for each top
        self.alist1 = [] # closest to pivot 1
        self.alist2 = [] # closest to pivot 2

        sets = [ [i] for i in range(natom) ]
        # print("SETS0", sets) # TODO
        # print(self.pivot1,self.pivot2)
        idx = self.bonds.index([min(self.pivot1,self.pivot2),
                                max(self.pivot1,self.pivot2)])

        # Remove the bond between the two tops to create disjoint set of
        # connected atoms. TODO
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

        # Remove the pivot atoms from the atomlists
        # self.alist1.remove(self.pivot1)
        # self.alist2.remove(self.pivot2)

        if len(self.alist1 + self.alist2) != self.masses.size:
            print("WARNING: ATOMS MISSING FROM ONE OR BOTH TOPS OF THIS ROTOR")

        # print("ALIST1", self.alist1, self.pivot1, len(self.alist1))
        # print("ALIST2", self.alist2, self.pivot2)
        return


    ######################################################################
    def read_scan_data(self, file_name):
        '''
        Read the atoms in the dihedral angle and the energies of the 1-D PES
        from the log file of the calculation that scans the dihedral angle.

        .. note:
            The energies are saved in J/mol (SI units).
        '''

        energy_offset = 0
        energies = []

        # Step through file one line at a time
        with open(file_name, 'r') as f:
            lines = f.readlines()

            # print(er.findall())

            parsing_rotor = False #Boolean for parsing PES energy data

            for line in lines:

                # Get pivot atoms
                if re.search('D \s+\d', line):
                    # print(line.split()) # TODO
                    splt = line.split()
                    self.dihedral = np.array([int(splt[1]), int(splt[2]),
                                              int(splt[3]), int(splt[4])]) - 1
                    self.pivot1 = self.dihedral[1]
                    self.pivot2 = self.dihedral[2]
                    self.scan_rate = float(splt[-1]) # In degrees

                # Make sure that none of the optimizations stopped prematurely
                if re.search('\s+Optimization stopped.', line):
                    # print(line.split())
                    print("Your rotational PES scan did not complete!")
                    print("We suggest that you increase the maximum number of optimiztion cycles or loosen your tolerances.")
                    print("Exiting...")
                    exit(0)

                # Start of scan data
                if re.match(' Summary of Optimized Potential Surface Scan', line):
                    if len(line.split()) == 6:
                        energy_offset = 0
                    elif len(line.split()) == 10:
                        energy_offset = float(line.split()[7])
                    else:
                        print("Something isn't right with your PES energies")
                        exit(0)
                    parsing_rotor = True

                # Energies
                if parsing_rotor:
                    if re.match('\s+Eigenvalues --', line):
                        # For G09
                        if len(line.split()) == 3:
                            for e in line.split()[-1].split('-')[1:]:
                                energies.append(- float(e))

                        # For G16
                        else:
                            for e in line.split()[2:]:
                                energies.append(float(e))

                # End of scan data
                if parsing_rotor and re.match(' GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad', line):
                    parsing_rotor = False
                    break

        # TODO Clean this up and just use self.* throughout function
        self.energy_offset = energy_offset * ha_to_j
        self.energies = np.array(energies)*ha_to_j - self.energy_offset # J


    ######################################################################
    def fit_potential(self, a=None, e=None, n_term=11):
        '''
        A is the matrix of fourier terms plus the constraint that
        :math:`\frac{dV}{d\theta} = 0` at :math:`\theta=0`. e are the energies.
        N_term is the number of Fourier terms to include in the expansion in
        units of J.
        '''

        if e == None: e = self.energies.copy()

        # E1 contains all the scan energies and the constraint.
        e1 = np.zeros((e.size+1)) #TODO
        e1[:-1] = e[:]
        e1[-1] = 0

        if a == None:
            a = np.zeros( (e.shape[0]+1, n_term) )
            th = np.linspace(0,2*np.pi,e.shape[0])

            for r in range(a.shape[0]):
                # Add constraint that dV/d\theta = 0 at \theta = 0
                if r == a.shape[0] -1:
                    for c in range(a.shape[1]):
                        if c == 0 or c%2 == 0:
                            a[r,c] = 0
                        else:
                            a[r,c] = int((c+1)/2)
                            # a[r,c] = 1
                            # print(a[r,c])
                    continue

                for c in range(a.shape[1]):
                    # First coefficient (no sinusoidal function)
                    if c == 0:
                        a[r,c] = 1
                    # a coeffs
                    elif c%2 == 0 and c > 0:
                        # print(int((c + 0)/2))
                        a[r,c] = np.cos(int((c + 0)/2) * th[r])
                    # b coeffs
                    else:
                        # print(int((c + 1)/2))
                        a[r,c] = np.sin(int((c + 1)/2) * th[r])

        # # print(a[1])
        # N = e.shape[0]
        # A = np.zeros((N+1,2*n_term), np.float64)
        # b = np.zeros(N+1, np.float64)
        # for i in range(N):
        #     phi = th[i]
        #     for m in range(n_term):
        #         A[i,m] = np.cos(m * phi)
        #         A[i,n_term+m] = np.sin(m * phi)
        #         b[i] = e[i]
        # # This row forces dV/dangle = 0 at angle = 0
        # for m in range(n_term):
        #     A[N,m+n_term] = 1
        # print(A[:,n_term])
        # print(a[-1,:])
        # print("A[N] ",A[N])
        # Add constraint to the e1


        pot_coeffs = reduce(np.dot, (la.inv(np.dot(a.T, a)), a.T, e1))
        # Store them in an easier to use manner
        self.v0_coeff = pot_coeffs[0] # V_0 coeff
        self.pot_coeffs = np.zeros((2,int((n_term-1)/2)))
        for i in range(self.pot_coeffs.shape[1]):
            self.pot_coeffs[0,i] = pot_coeffs[2*i+2] # a_k coeffs
            self.pot_coeffs[1,i] = pot_coeffs[2*i+1] # b_k coeffs


        # print(self.pot_coeffs[0])
        # print(self.pot_coeffs[1])
        # check1 = np.array([(k+1)*self.pot_coeffs[1,k] for k in range(0,int((n_term-1)/2))])
        # check2 = np.array([(k+1)**2*self.pot_coeffs[0,k] for k in range(0,int((n_term-1)/2))])
        # print(np.sum(check1))
        # print(np.sum(check2))
        # e = e[:-1] #TODO
        # print(e)

    ######################################################################
    def calc_reduced_moment(self):
        '''
        Calculate the reduced moment of inertia for the rotor.
        '''

        masses = self.masses
        geom = self.geom

        # Get the center of mass for each top
        m1 = 0.0
        m2 = 0.0
        cm1 = np.zeros((3,))
        cm2 = np.zeros((3,))

        for a1 in self.alist1:
            cm1 += masses[a1] * geom[a1,:]
            m1 += masses[a1]

        for a2 in self.alist2:
            cm2 += masses[a2] * geom[a2,:]
            m2 += masses[a2]

        cm1 /= m1
        cm2 /= m2

        # Moments of inertia
        ax_of_rot = (cm1-cm2)/la.norm(cm1-cm2)
        I1 = 0.0
        I2 = 0.0

        for a1 in self.alist1:
            r1 = (geom[a1,:] - cm1) - (np.dot((geom[a1,:] - cm1), ax_of_rot) *
                ax_of_rot)
            I1 += masses[a1] * la.norm(r1)**2

        for a2 in self.alist2:
            r2 = (geom[a2,:] - cm2) - (np.dot((geom[a2,:] - cm2), ax_of_rot) *
                ax_of_rot)
            I2 += masses[a2] * la.norm(r2)**2

        # Reduced moment of inertia
        self.i_red = 1.0/ (1.0/I1 + 1.0/I2)

        # Change units from g/mol*A^2 to kg*m^2
        self.i_red /= 1.0e23
        self.i_red /= N_avo


    ######################################################################
    def calc_e_levels(self, save_fig=False):
        '''
        Calculate the hindered rotor energies for the rotor object.

        Energy levels are in J
        '''

        m = 250
        ham_ir = np.zeros((m*2+1, m*2+1), dtype=np.complex_) # hindered rotor hamiltonian

        i_red = self.i_red
        coeffs = self.pot_coeffs
        v0 = self.v0_coeff
        # print(v0)
        n_fit = coeffs.shape[1] #int((pot_coeffs.size - 1)/2)
        B = h**2 / (8.0 * np.pi**2 * i_red)
        for n in range(-m, m+1):
            # (kg * m^2/s)^2 / (kg * m^2) [=] kg * m^2 / s^2 [=] J
            # print(n**2 * h**2 / (8.0 * np.pi**2 * i_red)) # TODO
            # print(v0)
            ham_ir[n+m,n+m] = n**2 *  B  + v0

            for k in range(1, n_fit+1):
                # rows > cols (below the diagonal)
                if n + m - k >= 0:
                    ham_ir[n+m,n+m-k] = 0.5 * coeffs[0,k-1] # a_k term
                    ham_ir[n+m,n+m-k] += 0.5j * coeffs[1,k-1] # b_k term
                # cols > rows (above the diagonal)
                if n + m + k < 2*m+1:
                    ham_ir[n+m,n+m+k] = 0.5 * coeffs[0,k-1] # a_k term
                    ham_ir[n+m,n+m+k] -= 0.5j * coeffs[1,k-1] # b_k term

        self.e_levels, _ = la.eigh(ham_ir) # TODO
        if save_fig:
            plt.figure()
            v_fit = self.calculate_potential()
            v_fit = np.roll(v_fit, 100)
            th = np.linspace(0,2*np.pi, num=v_fit.size)
            plt.plot(th, v_fit)

            for i in range(100):
                plt.plot(th, [self.e_levels[i]]*th.size)

            # plt.show()
            plt.savefig("e_levels_%.0f.png" % self.v_ho)
            # plt.close('all')

    ######################################################################
    ### Optional Functions ###############################################
    ######################################################################
    def calculate_thermo(self, t, harmonic=False):
        '''
        Calculate the contributions to the thermodynamic from the rigid rotor.

        Arguments:
            t: temperature in K

            harmonic: boolean flag for whether or not to calculate the
                properties harmonically

        Returns:
            q_ir: partition function

            h_ir: enthalpy in kcal/mol

            cp_ir: heat capacity in cal/(mol K)

            s_ir: entropy in cal/(mol K)


        TODO: Need to add degeneracy?
        '''

        if not harmonic:
            eps = self.e_levels.copy()
        else:
            mu = c_in_cm * self.v_ho
            eps = np.array([h*mu*(n+0.5) for n in range(10000)])

        zpe = eps[0] * N_avo/ 1.0e3/ kcal_to_kj
        # print(eps[0]/kb/t)
        eps -= eps[0] # Shift by ZPE
        eps *= N_avo # Convert to molar quantities (J/mol)

        bw = np.exp(-eps /(R * t)) # Boltzmann weights

        # Partition Function
        q_ir = np.sum(bw)
        q_ir /= self.sym

        # Enthalpy
        h_ir = np.sum(eps * bw) / np.sum(bw) /1.0e3 /kcal_to_kj # kcal/mol

        # Heat Capacity
        cp_ir = (np.sum(bw) * np.sum(eps*eps*bw) - np.sum(eps*bw)**2)/ \
            (np.sum(bw)**2 * R * t**2 * cal_to_j) # cal/(mol K)
        # print("T %i \t Cp_ir %f" %(t,cp_ir))

        # Entropy
        s_ir = (np.sum(eps*bw)/np.sum(bw) / t + \
            np.log(q_ir) * R )/ cal_to_j # cal/(mol K)

        return q_ir, h_ir, cp_ir, s_ir

    ######################################################################
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


    ######################################################################
    def write_potential_to_file(self, energy_offset, energies, file_name):
        pot_file_header = "# 1-D Scan of Total Energy\n#\n# %8s\t%16s\n"%('Step Number','Energy (Ha)')
        with open(file_name, 'w') as f:
            f.write(pot_file_header)
            i = 1;
            for energy in energies:
                f.write("%8i\t\t\t%16.8f\n"%(i,energy+energy_offset))
                i += 1


    ######################################################################
    def plot_rotational_pes(self, fitted=True, n_fit=100, show=False):
        '''
        Plot the 1-D Rotor PES.
        '''

        if fitted:
            th = np.linspace(0, 2*np.pi, num=n_fit)
            v_fit = self.calculate_potential(th=th)

        plt.figure()
        plt.plot(np.linspace(0,2*np.pi, num=len(self.energies)),
                 np.array(self.energies), marker='o', linestyle='--',
                 label="Calculated")
        plt.plot(th, v_fit, label="Fit")
        plt.xlabel('Angle / Radians')
        plt.ylabel('Energy / J')
        plt.legend()
        if show: plt.show()
        plt.savefig('rotor_%s.png' % self.label)
        # plt.close('all')



    ######################################################################
    def calculate_potential(self, th=None, n_fit=1000 ):
        '''
        Calculate the fitted potential using the pot_coeffs array.
        '''

        if th is None:
            th = np.linspace(0, 2*np.pi, num=n_fit)

        v_fit = np.zeros(th.shape[0])

        if self.pot_coeffs.shape[0] == 0: self.fit_potential()

        for i in range(th.shape[0]):
          v_fit[i] += self.v0_coeff * 1.0
          for j in range(self.pot_coeffs.shape[1]):
              v_fit[i] += self.pot_coeffs[0,j] * np.cos((j+1) * th[i])
              v_fit[i] += self.pot_coeffs[1,j] * np.sin((j+1) * th[i])

        return v_fit

    ######################################################################
    def print_properties(self, out_file):
        '''
        Write the properties of the rotors the specifies out_file.

        Arguments:

            out_file: open output file object to write the properties to
        '''

        # Write header
        header = 'Hindered Rotor with Harmonic Freq: %8.2f\n' % self.v_ho
        header += '-'*(len(header)-1) + '\n\n'
        out_file.write(header)

        # Print atoms in each top
        out_file.write('    Pivot atoms: %i %i\n' % (self.pivot1, self.pivot2))
        out_file.write('    Atoms in top 1: ')
        for atom in self.alist1: out_file.write('%i ' % atom)
        out_file.write('\n')
        out_file.write('    Atoms in top 2: ')
        for atom in self.alist2: out_file.write('%i ' % atom)
        out_file.write('\n\n')

        # Fourier Series coefficients
        out_file.write('    Fourier series coefficients (J):\n')
        out_file.write('    V_0 = %8.3e\n' % self.v0_coeff)
        out_file.write('    %16s %16s\n' % ('Cosine Coeffs.', 'Sine Coeffs.'))
        for i in range(self.pot_coeffs[0].size):
            out_file.write('    %16.3e %16.3e\n' % (self.pot_coeffs[0,i],
                                                self.pot_coeffs[1,i]))
        out_file.write('\n\n')

        # Reduced moment of inertia
        out_file.write('    Reduced moment of inertia (kg*m^2): %8.3e' %
                       self.i_red)
        out_file.write('\n\n')


    ############################################################################
    ### Deprecated #############################################################
    ############################################################################

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

################################################################################

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
    k /= la.norm(k)

    # Angle we need to rotate
    th = np.arccos( np.dot(v1,v2)/(la.norm(v1)*la.norm(v2)) )

    # Euler/Rodrigues params
    # See https://en.wikipedia.org/wiki/Euler%E2%80%93Rodrigues_formula
    a = np.cos(th/2.)
    b = k[0] * np.sin(th/2.)
    c = k[1] * np.sin(th/2.)
    d = k[2] * np.sin(th/2.)

    print("CHECK %f = 1" % (a**2 + b**2 + c**2 + d**2))

    r = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                  [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                  [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

    return r
