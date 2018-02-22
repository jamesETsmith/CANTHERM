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

import re #, os, sys
import numpy as np
import matplotlib.pyplot as plt

#import pandas as pd
#import seaborn as sns

def read_scan_data(file_name):
    energy_offset = 0
    energies = []

    with open(file_name, 'r') as f:
        lines = f.readlines()
        parsing_rotor = False


        for line in lines:
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

    return energy_offset, energies

def write_potential_to_file(energy_offset, energies, file_name):
    pot_file_header = "# 1-D Scan of Total Energy\n#\n# %8s\t%16s\n"%('Step Number','Energy (Ha)')
    with open(file_name, 'w') as f:
        f.write(pot_file_header)
        i = 1;
        for energy in energies:
            f.write("%8i\t\t\t%16.8f\n"%(i,energy+energy_offset))
            i += 1

    return


if __name__ == "__main__":

    energy_offset1, energies1 = read_scan_data('examples/h-transf/oxy_rotor/oxy_rotor.log')
    energy_offset2, energies2 = read_scan_data('examples/h-transf/methyl_rotor/methyl_rotor.log')
    energy_offset3, energies3 = read_scan_data('examples/h-transf/hydroxyl_rotor/hydroxyl_rotor.log')

    write_potential_to_file(energy_offset1, energies1, 'examples/h-transf/1.txt')
    write_potential_to_file(energy_offset2, energies2, 'examples/h-transf/2.txt')
    write_potential_to_file(energy_offset3, energies3, 'examples/h-transf/3.txt')

    plt.figure()
    plt.plot(energies1, label='oxy_rotor')
    plt.plot(energies2, label='methyl_rotor')
    plt.plot(energies3, label='hydroxy_rotor')
    plt.legend()
    plt.show()
