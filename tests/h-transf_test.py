#!/usr/bin/env python
'''
Translations
Q 2.555e+07
Cv cal/(mol K) 4.968016
H kcal/mol 1.481214
S cal/(mol K) 38.862303

Rotations
Q 8.022e+04 # Use 0.884719D+05 from Gaussian
Cv cal/(mol K) 2.980810
H kcal/mol 0.888728
S cal/(mol K) 25.421505 # Use 25.616 from Gaussian

Vibrations
Q 1.110e+01
Cv cal/(mol K) 15.883897
H kcal/mol 2.170924
S cal/(mol K) 12.065042
'''

'''
2.555e+07        1.480          4.965          38.836
'''

import unittest
import os
import cantherm.cantherm as cantherm

sig_fig = 4

os.chdir('../examples/h-transf/')
data = cantherm.CanTherm(input_filename = 'input2', verbose=0)
data.run()
gs = data.MoleculeList[0]

class KnowValues(unittest.TestCase):
    def test1_trans(self):
        print('\nTesting translational thermodynamics...')
        q_tr, h_tr, cp_tr, s_tr = gs.calculate_thermo(False, [298.15],
                                                           mode_type='trans')
        self.assertAlmostEqual(q_tr[0]/1e7, 2.5552579, places=(sig_fig-1) )
        self.assertAlmostEqual(h_tr[0], 1.4802253, places=(sig_fig-1) )
        self.assertAlmostEqual(cp_tr[0], 4.9647000, places=(sig_fig-1) )
        self.assertAlmostEqual(s_tr[0]/10, 3.8836363, places=(sig_fig-1) )

    def test2_rot(self):
        print('\nTesting rotational thermodynamics...')
        q_rot, h_rot, cp_rot, s_rot = gs.calculate_thermo(False, [298.15],
                                                           mode_type='rot')
        self.assertAlmostEqual(q_rot[0]/1e4, 8.8472248, places=(sig_fig-1))
        self.assertAlmostEqual(h_rot[0]*10, 8.88135182, places=(sig_fig-1) )
        self.assertAlmostEqual(cp_rot[0], 2.9788200, places=(sig_fig-1) )
        self.assertAlmostEqual(s_rot[0]/10, 2.5598875, places=(sig_fig-1) )

    def test3_vib(self):
        print('\nTesting vibrational thermodynamics...')
        q_vib, h_vib, cp_vib, s_vib = gs.calculate_thermo(False, [298.15],
                                                           mode_type='vib')
        self.assertAlmostEqual(q_vib[0]/10, 1.1103774, places=(sig_fig-1) )
        self.assertAlmostEqual(h_vib[0], 2.1709434, places=(sig_fig-1) )
        self.assertAlmostEqual(cp_vib[0]/10, 1.5873404, places=(sig_fig-1) )
        self.assertAlmostEqual(s_vib[0]/10, 1.2057099, places=(sig_fig-1) )

    def test4_kin(self):
        print('\nTesting kinetics...')
        rx = data.rx

        # TST Rate Tests
        self.assertAlmostEqual(rx.rates[0]*10, 1.0079000, places=(sig_fig-1) )
        self.assertAlmostEqual(rx.rates[1]*10, 1.2223926, places=(sig_fig-1) )
        self.assertAlmostEqual(rx.rates[2]/1e2, 2.7727886, places=(sig_fig-1) )
        self.assertAlmostEqual(rx.rates[3]/1e4, 2.7434532, places=(sig_fig-1) )
        self.assertAlmostEqual(rx.rates[4]/1e5, 5.8109805, places=(sig_fig-1) )
        self.assertAlmostEqual(rx.rates[5]/1e7, 2.6668462, places=(sig_fig-1) )
        self.assertAlmostEqual(rx.rates[6]/1e8, 2.7148873, places=(sig_fig-1) )
        self.assertAlmostEqual(rx.rates[7]/1e9, 6.4290564, places=(sig_fig-1) )

        # Test Fitted Arrhenius Params
        self.assertAlmostEqual(rx.arrhenius_coeff/1e12, 2.678,
                               places=(sig_fig-1) )
        self.assertAlmostEqual(rx.arrhenius_exp, 0.125630, places=(sig_fig-1) )
        self.assertAlmostEqual(rx.activation_energy, 18.203,
                               places=(sig_fig-1) )

if __name__ == "__main__":
    print("Full Tests for H-Transf")
    unittest.main()
