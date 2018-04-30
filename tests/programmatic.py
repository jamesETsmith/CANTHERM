#!/usr/bin/env python

import unittest
import os
# import cantherm.main as main
import cantherm.molecule as molecule

sig_fig = 4
gs = molecule.Molecule('../examples/h-transf/gs.log', False, 0.99, 0, './', txt_input=False)


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

    def test4_total(self):
        print('\nTesting total thermodynamics...')
        q_tot, h_tot, cp_tot, s_tot = gs.calculate_all_thermo([298.15], False)
        self.assertAlmostEqual(q_tot[0]/1e13, 2.510, places=(sig_fig-1) )
        self.assertAlmostEqual(h_tot[0], 4.539, places=(sig_fig-1) )
        self.assertAlmostEqual(cp_tot[0]/10, 2.3817, places=(sig_fig-1) )
        self.assertAlmostEqual(s_tot[0]/10, 7.6492, places=(sig_fig-1) )


if __name__ == "__main__":
    print("Programmatic Tests for H-Transf")
    unittest.main()
