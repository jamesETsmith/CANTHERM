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

import math, re
import numpy as np
import matplotlib.pyplot as plt
from cantherm.constants import *
from functools import reduce
from cantherm import readGeomFc
readEnergy = readGeomFc.readEnergy

la = np.linalg

class Reaction:
    '''
    Attributes:

        reactants ():

        ts ():

        temp ():

        calc_type ():

        reac_type ():

        products ():

        tunneling ():

        scale ():

        rates ():

        arrhenius_coeff ():

        activation_energy ():

    '''

    def __init__(self, reactants, ts, temp, calc_type='CBS-QB3',
                 reac_type='Unimol', products=[], tunneling=None, scale=0.99):
        self.reactants = reactants
        self.ts = ts
        self.temp = temp
        self.calc_type = calc_type
        self.reac_type = reac_type
        self.products = products
        self.tunneling = tunneling
        self.scale = scale
        self.rates = None
        self.tunneling_coeff = None
        self.q_ratio = None
        self.arrhenius_coeff = None
        self.activation_energy = None
        return

    ############################################################################

    def calc_TST_rates(self):
        '''
            Calculates the transition state theory rate constants for the
            reaction at the given temperatures.
        '''
        self.rates = [0]*len(self.temp)
        self.tunneling_coeff = [0]*len(self.temp)
        self.q_ratio = [0]*len(self.temp)

        for i in range(len(self.temp)):
            t = self.temp[i]
            if self.reac_type == 'Unimol':
                Q_react = self.reactants.calculate_Q(t)
                Q_TS = self.ts.calculate_Q(t)

            self.q_ratio[i] = (kb * t/ h) * (Q_TS/ Q_react)
            self.rates[i] = (kb * t/ h) * (Q_TS/ Q_react)
            self.rates[i] *= math.exp(-(self.ts.Energy - \
                                        self.reactants.Energy) * ha_to_kcal \
                                        / R_kcal / t)


            if self.tunneling == "Wigner":
                kappa = wigner_correction(t, self.ts.imagFreq, self.scale)
                self.rates[i] *= kappa
                self.tunneling_coeff[i] = kappa

            # print("Test rate constant for T=%i \t %e" % (T,k_test*kappa[i]))


    ############################################################################

    def fit_arrhenius(self):
        '''
        Fit the rates calculated using TST to an Arrhenius plot and return the
        A coefficient and the activation energy.

        .. note::
            This is a modified version :math:`k = A \frac{T}{T_0}^n \exp{-\frac{E_a}{R T}}`

        This boils down to solving the least squares problem Cx=b where the
        first column of C is 1, the second column is the log of the current
        temperature divided by the minimum temperature, and the third column is
        -1/RT for each T in self.temp. Where x[0] = ln(A), x[1] = E_a, and
        b[i]=ln(rate[i]). This can be solved using the normal equations.
        '''

        # In case the rates haven't been calculated yet.
        if self.rates == None: self.calc_TST_rates()

        # Construct C
        c = np.ones( (len(self.temp),3) )
        c[:,1] = np.log(np.array(self.temp)/1000)
        c[:,-1] = 1./np.array(self.temp)
        b = np.log(self.rates)

        x = reduce(np.dot, (la.inv(np.dot(c.T,c)), c.T, b ))

        self.arrhenius_coeff = np.exp(x[0])
        self.arrhenius_exp = x[1]
        self.activation_energy = - x[-1] * R_kcal


    ############################################################################

    def print_arrhenius(self):
        '''
            Print out the fitted Arrhenius plot.

            :math:`\log{k} = \log{A} + n \log{T/1000} - \frac{E_a}{R T}`
        '''

        # In case the fitting hasn't taken place yet.
        if self.arrhenius_coeff == None:
            self.fit_arrhenius()

        t_0 = 1000 #np.min(self.temp) # Minimum temp in user specified range
        t_inv = 1.0/np.linspace(self.temp[0], self.temp[-1], 1000)

        fit_k = np.log(1/(t_inv * t_0)) * self.arrhenius_exp
        fit_k += np.log(self.arrhenius_coeff)
        fit_k -= self.activation_energy*t_inv/R_kcal

        plt.figure()
        plt.scatter(1./np.array(self.temp), np.log(self.rates),label="TST Data")
        plt.plot(t_inv, fit_k, 'r', label="Fitted Data")
        plt.legend()
        plt.show()


    ############################################################################

    def print_properties(self, out_file):
        '''
        Write the reaction rates and fitted arrhenius coefficients to the output
        file.
        '''

        # Write TST header
        kin_header = '%9s   %14s   %14s   %15s\n' % ('Temp. (K)',
                                                    'Q Ratio (s^-1)',
                                                    'k(TST) (s^-1)',
                                                    'k(TST+T) (s^-1)')
        kin_header += '-'*9 + '   ' + '-'*14 + '   ' + '-'*14 + '   ' + '-'*15 + '\n\n'
        out_file.write(kin_header)

        for i in range(len(self.temp)):
            out_file.write("%9i   %14.3e   %14.3e   %15.3e\n"%(self.temp[i],
                self.q_ratio[i],
                self.rates[i]/self.tunneling_coeff[i],
                self.rates[i]))
        out_file.write('\n')

        # Arrhenius Data
        arr_header = 'Fitted Arrhenius Data\n'
        arr_header += '-'*(len(arr_header) -1) + '\n\n'
        arr_eqn = '\log{k} = \log{A} + n \log{T/1000} - frac{E_a}{R T}\n\n'
        out_file.write(arr_header + arr_eqn)
        out_file.write('A = %.3e s^1\n' % self.arrhenius_coeff)
        out_file.write('n = %.6f \n' % self.arrhenius_exp)
        out_file.write('E_a = %.3f kcal/mol\n' % self.activation_energy)



################################################################################

class RxSystem:

    def __init__(self, outputs, temps, name):
        '''
        Arguments:

            outputs ([string]): List of reaction output files.

            temps ([float]): List of temperatures.

            name (string): Name of the rxn system
        '''

        self.outputs = outputs
        self.temps = temps
        self.name = name
        self.a_coeffs = []
        self.a_exp = []
        self.ea = []
        self.q_ratios = []
        self.tst_rates = []
        self.tstt_rates = []
        self.energy_files = []
        self.barriers = []


    ############################################################################

    def get_kinetics_from_output(self):

        n_temps = len(self.temps)
        for filename in self.outputs:
            with open(filename, 'r') as f:
                reading_kinetics = False
                lines = f.readlines()
                rates_tstt = []
                rates_tst = []
                q_ratio = []
                e_files = []

                ctr = 0
                for i in range(len(lines)):
                    if re.search('\s+Kinetic Properties*', lines[i]):
                        reading_kinetics = True
                        # i += 7
                        continue

                    if reading_kinetics and ctr >= n_temps:
                        reading_kinetics = False
                        continue

                    if reading_kinetics:
                        if len(lines[i].split()) == 4 and lines[i].split()[0] != '---------':
                            # print(lines[i].split())
                            q_ratio.append(float(lines[i].split()[1]))
                            rates_tst.append(float(lines[i].split()[2]))
                            rates_tstt.append(float(lines[i].split()[3]))
                            ctr += 1

                    if len(lines[i].split(':')) > 0 and \
                        lines[i].split(':')[0] == 'Energy file':
                        # Add the file name and the energy type
                        e_files.append([lines[i].split(':')[1].strip(),
                                        lines[i+1].split()[-1]])

                self.energy_files.append(e_files)
                # print(e_file)
                # print(rates)
                self.q_ratios.append(q_ratio)
                self.tst_rates.append(rates_tst)
                self.tstt_rates.append(rates_tstt)



    ############################################################################

    def plot_tst_rates(self, log_plot=True):
        plt.figure()
        i = 0

        for output in self.outputs:
            if log_plot:
                rates = np.log(np.array(self.tst_rates[i]))
            else:
                rates = self.tst_rates[i]
            plt.plot(self.temps, rates, 'o', ls='-', label=output)

            i += 1

        plt.legend()
        plt.xlabel("Temperature / K")
        if log_plot:
            plt.ylabel("ln(k)")
        else:
            plt.ylabel("k / $s^{-1}$")
        plt.savefig(self.name+".png")

    ############################################################################

    def tabulate_tst_rates(self, comp_data=[], comp_labels=[]):

        plt.figure()

        # plt.axis('tight')
        plt.axis('off')
        plt.grid('off')

        # Columns are the Temperatures and Rows are the specific reactions
        cellText = []
        for i in range(len(self.tst_rates)):
            cellText.append(["%.3e" % x for x in self.tstt_rates[i]])

        # If applicatble append comparison data/labels
        for i in range(len(comp_data)):
            row = ["%.3e" % x for x in comp_data[i]]
            row += ["-"]*(len(self.temps)-len(comp_data))
            cellText.append(row)

        # Row labels
        rowLabels = self.outputs + comp_labels
        if len(cellText) != len(rowLabels):
            rowLabels += ['Comparison']*(len(cellText)-len(rowLabels))

        # Column labels
        colLabels = ['%d (K)' % t for t in self.temps]

        # Add a table at the bottom of the axes
        the_table = plt.table(cellText=cellText, rowLabels=rowLabels,
                              colLabels=colLabels, loc='best')

        plt.savefig(self.name+"_table"+".pdf", bbox_inches="tight")

################################################################################

    def tabulate_kinetics(self, labels, comp_data):

        plt.figure()

        # plt.axis('tight')
        plt.axis('off')
        plt.grid('off')

        # Columns are the Temperatures and Rows are the specific reactions
        cellText = []
        for i in range(len(self.tst_rates)):
            row = ['%.3e'%self.tst_rates[i][0],
                    '%.3e'%self.tstt_rates[i][0],
                    '%.3e'%comp_data]
            if i == 0:
                row[-1] = '-'
            cellText.append(row)


        # Row labels
        rowLabels = labels

        # Column labels
        colLabels = ['TST (s^-1)', 'TST+T (s^-1)', 'V&P (s^-1)']

        # Add a table at the bottom of the axes
        the_table = plt.table(cellText=cellText, rowLabels=rowLabels,
                              colLabels=colLabels, loc='best')

        plt.savefig(self.name+"_table"+".pdf", bbox_inches="tight")

################################################################################

def get_barriers(energy_files, outputs, plot=False, plot_name='', comp_data=[],
                 comp_labels=[]):
    '''
    Collect the energy barrier data from the various final and possibly
    other intermediate eneriges, e.g. MP2 energies in a CCSD(T) calculation.

    Arguments:

        energy_files ([[[string]]]): First index is reaction, second is gs or ts,
            third is energy (0) and method (1).
    '''

    barriers = [[],[],[]]

    i = 0
    for e_files in energy_files:
        if e_files[0][1] != e_files[1][1]:
            print("ERROR: Energy methods for GS and TS don't match.")
            print(e_files[0][1], e_files[1][1])
            print("Exiting...")
            exit(0)

        gs_e = readEnergy(e_files[0][0], e_files[0][1])
        ts_e = readEnergy(e_files[1][0], e_files[1][1])
        # print(gs_e, ts_e, e_files[0][1], e_files[1][1])

        # Barrier
        dE = ts_e - gs_e
        # Add barrier along with method to barriers attribute.
        barriers[0].append(dE)
        barriers[1].append(e_files[0][1])
        barriers[2].append( outputs[i])

        # Get "intermediate" energies
        if e_files[0][1] == 'cbsqb3':
            gs_ub3lyp_e = readEnergy(e_files[0][0], 'ub3lyp')
            ts_ub3lyp_e = readEnergy(e_files[1][0], 'ub3lyp')

            dE_ub3lyp = ts_ub3lyp_e - gs_ub3lyp_e
            barriers[0].append(dE_ub3lyp)
            barriers[1].append('ub3lyp')
            barriers[2].append( outputs[i])

        elif e_files[0][1] == 'DF-LUCCSD(T)-F12':
            gs_lrmp2_e = readEnergy(e_files[0][0], 'RHF-LRMP2')
            ts_lrmp2_e = readEnergy(e_files[1][0], 'RHF-LRMP2')

            dE_lrmp2 = ts_lrmp2_e - gs_lrmp2_e
            barriers[0].append(dE_lrmp2)
            barriers[1].append('RHF-LRMP2')
            barriers[2].append( outputs[i])

            gs_ccsd_e = readEnergy(e_files[0][0], 'LUCCSD-F12a')
            ts_ccsd_e = readEnergy(e_files[1][0], 'LUCCSD-F12a')

            dE_ccsd = ts_ccsd_e - gs_ccsd_e
            barriers[0].append(dE_ccsd)
            barriers[1].append('LUCCSD-F12a')
            barriers[2].append( outputs[i])

        elif e_files[0][1] == 'UCCSD(T)-F12':
            gs_rmp2_e = readEnergy(e_files[0][0], 'RHF-RMP2')
            ts_rmp2_e = readEnergy(e_files[1][0], 'RHF-RMP2')

            dE_rmp2 = ts_rmp2_e - gs_rmp2_e
            barriers[0].append(dE_rmp2)
            barriers[1].append('RHF-RMP2')
            barriers[2].append( outputs[i])

            gs_ccsd_e = readEnergy(e_files[0][0], 'RHF-UCCSD-F12a')
            ts_ccsd_e = readEnergy(e_files[1][0], 'RHF-UCCSD-F12a')

            dE_ccsd = ts_ccsd_e - gs_ccsd_e
            barriers[0].append(dE_ccsd)
            barriers[1].append('RHF-UCCSD-F12a')
            barriers[2].append( outputs[i] )

        i += 1

    # print(barriers)
    # Plot the barrier heigts
    if plot:
        all_bars = barriers[0] + comp_data
        ind = np.arange(len(all_bars))
        wid = 0.35
        colors = ['blue']*(len(barriers[0])) + ['red']*len(comp_data)

        plt.figure()
        plt.bar(ind, all_bars, wid, align='edge', color=colors)

        # Label the height of the bars
        for i in range(len(all_bars)):
            height = all_bars[i]
            plt.text(i + wid/4., 1.05*height, '%.5f' % height, ha='center')

        # # Label the height of the bars for comparison data
        # for i in range(len(comp_data)):
        #     height = comp_data[i]
        #     plt.text(i+len(barriers[0]) + wid/4., 1.05*height,
        #              '%.5f' % height, ha='center')

        # Titles and ticks
        # If no labels are present just say comparison
        if len(comp_labels) != len(comp_data):
            comp_labels += ['Comparison']*(len(comp_data)-len(comp_labels))

        all_labels = barriers[1] + comp_labels

        plt.title('Barrier Heights as a Function of Method for %s' % plot_name)
        plt.xticks(ind+wid/2., all_labels, rotation=45)
        plt.ylabel('Energy Barrier / Ha')

        # Formatting
        plt.ylim( min(0,min(all_bars)*1.25), max(all_bars)*1.25 )
        plt.xlim( ind[0]-0.5, ind[-1]+1 )
        plt.tight_layout()
        plt.axhline(0, color='black')
        plt.savefig(plot_name+ '_barrier_heights'+'.png')
        plt.close()

    return(barriers)


################################################################################

def wigner_correction(t, freq, scale):
    # see doi:10.1103/PhysRev.40.749 and doi:10.1039/TF9595500001
    return (1.0 + 1.0 / 24.0 * (h * abs(freq) * scale * c_in_cm / (t * kb) )**2)
