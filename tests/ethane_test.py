'''
Translations
Q 6.474e+06
Cv cal/(mol K) 4.968016
H kcal/mol 1.481214
S cal/(mol K) 36.133899

Rotations
Q 8.091e+02
Cv cal/(mol K) 2.980810
H kcal/mol 0.888728
S cal/(mol K) 16.287065

Vibrations
Q 1.092e+00
Cv cal/(mol K) 3.084699
H kcal/mol 0.211794
S cal/(mol K) 0.886107

All Hindered Rotors
Q 1.336e+00
Cv cal/(mol K) 2.083101
H kcal/mol 0.339318
S cal/(mol K) 1.714244
'''


import unittest
import os
import cantherm.cantherm as cantherm

sig_fig = 4

os.chdir('../examples/ethane/')
data = cantherm.CanTherm(input_filename = 'input', verbose=0)
data.run()
gs = data.MoleculeList[0]

class KnowValues(unittest.TestCase):
    def test1_trans(self):
        print('\nTesting translational thermodynamics...')
        q_tr, h_tr, cp_tr, s_tr = gs.calculate_thermo(False, [298.15],
                                                           mode_type='trans')
        # print("%.7e %.7f %.7f %.7f"%(q_tr[0],h_tr[0],cp_tr[0],s_tr[0]))
        self.assertAlmostEqual(q_tr[0]/1e6, 6.4737593, places=(sig_fig-1) )
        self.assertAlmostEqual(h_tr[0], 1.4802253, places=(sig_fig-1) )
        self.assertAlmostEqual(cp_tr[0], 4.9647000, places=(sig_fig-1) )
        self.assertAlmostEqual(s_tr[0]/10, 3.6109787, places=(sig_fig-1) )

    def test2_rot(self):
        print('\nTesting rotational thermodynamics...')
        q_rot, h_rot, cp_rot, s_rot = gs.calculate_thermo(False, [298.15],
                                                           mode_type='rot')
        # print("%.7e %.7f %.7f %.7f"%(q_rot[0],h_rot[0],cp_rot[0],s_rot[0]))
        self.assertAlmostEqual(q_rot[0]/1e2, 8.0913485, places=(sig_fig-1))
        self.assertAlmostEqual(h_rot[0]*10, 8.881352, places=(sig_fig-1) )
        self.assertAlmostEqual(cp_rot[0], 2.9788200, places=(sig_fig-1) )
        self.assertAlmostEqual(s_rot[0]/10, 1.6276204, places=(sig_fig-1) )

    def test3_vib(self):
        print('\nTesting vibrational thermodynamics...')
        q_vib, h_vib, cp_vib, s_vib = gs.calculate_thermo(False, [298.15],
                                                           mode_type='vib')
        # print("%.7e %.7f %.7f %.7f"%(q_vib[0],h_vib[0],cp_vib[0],s_vib[0]))
        self.assertAlmostEqual(q_vib[0], 1.0601616, places=(sig_fig-1) )
        self.assertAlmostEqual(h_vib[0]*10, 1.581894, places=(sig_fig-1) )
        self.assertAlmostEqual(cp_vib[0], 2.5463732, places=(sig_fig-1) )
        self.assertAlmostEqual(s_vib[0]*10, 6.462335, places=(sig_fig-1) )

    def test4_ir(self):
        print('\nTesting internal rotational thermodynamics...')
        q_ir, h_ir, cp_ir, s_ir = gs.calculate_thermo(False, [298.15],
                                                           mode_type='int rot')
        # print("%.7e %.7f %.7f %.7f"%(q_ir[0],h_ir[0],cp_ir[0],s_ir[0]))
        self.assertAlmostEqual(q_ir[0], 1.3919058, places=(sig_fig-1) )
        self.assertAlmostEqual(h_ir[0]*10, 3.259944, places=(sig_fig-1) )
        self.assertAlmostEqual(cp_ir[0], 2.0110415, places=(sig_fig-1) )
        self.assertAlmostEqual(s_ir[0], 1.7505080, places=(sig_fig-1) )

if __name__ == "__main__":
    unittest.main()
