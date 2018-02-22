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

import math
import numpy as np
import matplotlib.pyplot as plt
from constants import *
from functools import reduce

la = np.linalg

class Reaction:
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
        self.arrhenius_coeff = None
        self.activation_energy = None
        return

    def calc_TST_rates(self):
        '''
            Calculates the transition state theory rate constants for the
            reaction at the given temperatures.
        '''
        self.rates = [0]*len(self.temp)
        for i in range(len(self.temp)):
            t = self.temp[i]
            if self.reac_type == 'Unimol':
                Q_react = self.reactants.calculate_Q(t)
                Q_TS = self.ts.calculate_Q(t)

            self.rates[i] = (kb * t/ h) * (Q_TS/ Q_react)
            self.rates[i] *= math.exp(-(self.ts.Energy - \
                                        self.reactants.Energy) * ha_to_kcal \
                                        * 1.e3 / R_kcal / t)

            if self.tunneling == "Wigner":
                kappa = wigner_correction(t, self.ts.imagFreq, self.scale)
                self.rates[i] *= kappa

            # print("Test rate constant for T=%i \t %e" % (T,k_test*kappa[i]))
        return

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
        c[:,1] = np.log(np.array(self.temp)/np.min(self.temp))
        c[:,-1] = 1./np.array(self.temp)
        b = np.log(self.rates)

        x = reduce(np.dot, (la.inv(np.dot(c.T,c)), c.T, b ))

        self.arrhenius_coeff = np.exp(x[0])
        self.arrhenius_exp = x[1]
        self.activation_energy = - x[-1] * R
        return

    def print_arrhenius(self):
        '''
            Print out the fitted Arrhenius plot.

            :math:`\log{k} = \log{A} + n \log{T/T_0} - \frac{E_a}{R T}`
        '''

        # In case the fitting hasn't taken place yet.
        if self.arrhenius_coeff == None:
            self.fit_arrhenius()

        t_0 = np.min(self.temp) # Minimum temp in user specified range
        t_inv = 1.0/np.linspace(self.temp[0], self.temp[-1], 1000)

        fit_k = np.log(1/(t_inv * t_0)) * self.arrhenius_exp
        fit_k += np.log(self.arrhenius_coeff)
        fit_k -= self.activation_energy*t_inv/R

        plt.figure()
        plt.scatter(1./np.array(self.temp), np.log(self.rates),label="TST Data")
        plt.plot(t_inv, fit_k, 'r', label="Fitted Data")
        plt.legend()
        plt.show()



################################################################################

def wigner_correction(t, freq, scale):
    # see doi:10.1103/PhysRev.40.749 and doi:10.1039/TF9595500001
    return (1.0 + 1.0 / 24.0 * (h * abs(freq) * scale * c_in_cm / (t * kb) )**2)
