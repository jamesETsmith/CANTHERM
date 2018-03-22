#!/usr/bin/env python
'''
Translations
Q 5.493e+07
Cv cal/(mol K) 4.968016
H kcal/mol 1.481214
S cal/(mol K) 40.383014

Rotations
Q 4.275e+05
Cv cal/(mol K) 2.980810
H kcal/mol 0.888728
S cal/(mol K) 28.746466


'''


import unittest
import os
import cantherm.cantherm as cantherm

sig_fig = 4

os.chdir('../examples/pvc/')
data = cantherm.CanTherm(input_filename = 'input', verbose=0)
data.run()
gs = data.MoleculeList[0]

class KnowValues(unittest.TestCase):
    def test1_trans(self):
        print('\nTesting translational thermodynamics...')
        q_tr, h_tr, cp_tr, s_tr = gs.calculate_thermo(False, [298.15],
                                                           mode_type='trans')
        self.assertAlmostEqual(q_tr[0]/1e7, 5.4926419, places=(sig_fig-1) )
        self.assertAlmostEqual(h_tr[0], 1.4802253, places=(sig_fig-1) )
        self.assertAlmostEqual(cp_tr[0], 4.9647000, places=(sig_fig-1) )
        self.assertAlmostEqual(s_tr[0]/10, 4.0356070, places=(sig_fig-1) )

    def test2_rot(self):
        print('\nTesting rotational thermodynamics...')
        q_rot, h_rot, cp_rot, s_rot = gs.calculate_thermo(False, [298.15],
                                                           mode_type='rot')
        self.assertAlmostEqual(q_rot[0]/1e5, 4.2752503, places=(sig_fig-1))
        self.assertAlmostEqual(h_rot[0]*10, 8.88135182, places=(sig_fig-1) )
        self.assertAlmostEqual(cp_rot[0], 2.9788200, places=(sig_fig-1) )
        self.assertAlmostEqual(s_rot[0]/10, 2.8727279, places=(sig_fig-1) )

    def test3_vib(self):
        print('\nTesting vibrational thermodynamics...')
        q_vib, h_vib, cp_vib, s_vib = gs.calculate_thermo(False, [298.15],
                                                           mode_type='vib')
        self.assertAlmostEqual(q_vib[0], 7.2996385, places=(sig_fig-1) )
        self.assertAlmostEqual(h_vib[0], 2.0972722, places=(sig_fig-1) )
        self.assertAlmostEqual(cp_vib[0]/10, 1.6416168, places=(sig_fig-1) )
        self.assertAlmostEqual(s_vib[0]/10, 1.0977172, places=(sig_fig-1) )


    def test4_ir(self):
        print('\nTesting internal rotational thermodynamics...')
        q_ir, h_ir, cp_ir, s_ir = gs.calculate_thermo(False, [298.15],
                                                           mode_type='int rot')
        self.assertAlmostEqual(q_ir[0]/10, 3.1111082, places=(sig_fig-1) )
        self.assertAlmostEqual(h_ir[0], 1.3827954, places=(sig_fig-1) )
        self.assertAlmostEqual(cp_ir[0], 5.0372122, places=(sig_fig-1) )
        self.assertAlmostEqual(s_ir[0]/10, 1.1469068, places=(sig_fig-1) )


if __name__ == "__main__":
    print("Full Tests for PVC")
    unittest.main()
