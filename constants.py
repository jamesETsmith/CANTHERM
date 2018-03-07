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
### Constants ###
R_kcal = 1.985 # in kcal mol^-1 K^-1 * 1000
R = 8.314472 # J mol^-1 K^-1
R_kJ = R / 1000
kb = 1.3806504e-23 #J/K
kb_in_cm = 0.695039 #cm^-1/K
h = 6.62606896e-34
hbar = 1.054571800e-34 #J s
amu = 1.660538921e-27
c_in_cm = 29979245800
N_avo = 6.0221409e23


### Energy Conversions ###
ha_to_kcal = 627.5095 # kcal/mol/Ha
ha_to_kj = 2625.50 # kJ/mol/Ha
kcal_to_kj = 4.184 # kJ/kcal
cal_to_j = 4.184 # J/cal
ha_to_j = 4.359744650e-18 # J/ha
