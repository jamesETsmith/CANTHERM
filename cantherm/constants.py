#!/usr/bin/env python
LICENSE = """
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

### Constants ### 4.184
R = 8.314462618  # J mol^-1 K^-1 From NIST
R_cal = R / 4.184  # cal/(mol K)
R_kcal = R_cal / 1.0e3  # in kcal mol^-1 K^-1
R_kJ = R / 1000
kb = 1.38064852e-23  # J K^-1
kb_in_cm = 0.695039  # cm^-1/K
h = 6.62607004e-34  # J s
hbar = 1.0545718e-34  # J s
amu = 1.660538921e-27
c_in_cm = 29979245800
N_avo = 6.0221409e23


### Energy Conversions ###
ha_to_kcal = 627.5095  # kcal/mol/Ha
ha_to_kj_mol = 2625.50  # kJ/mol/Ha
kcal_to_kj = 4.184  # kJ/kcal
kj_to_kcal = 1.0 / kcal_to_kj
cal_to_j = 4.184  # J/cal
j_to_cal = 1 / cal_to_j  # cal/J
ha_to_j = 4.359744650e-18  # J/Ha
# ha_to_kj = ha_to
ha_to_ev = 27.211386  # eV/Ha
kcal_to_ev = 0.043363  # eV mol/kcal
cal_to_ev = 0.000043363  # eV mol/cal

### Other Conversions ###
b0_to_m = 5.29177e-11  # m/bohr
b0_to_angst = 5.29177e-1  # angstroms/bohr
