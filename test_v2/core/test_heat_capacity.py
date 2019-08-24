import numpy as np
import pytest

from cantherm.constants import N_avo
from cantherm.core import cp_tr, cp_rot, cp_vib

npt = np.testing

ethane_freqs = np.array(
    [
        827.959,
        827.985,
        997.484,
        1219.445,
        1219.455,
        1410.060,
        1425.772,
        1505.088,
        1505.109,
        1507.738,
        1507.748,
        3025.133,
        3025.796,
        3071.286,
        3071.296,
        3096.595,
        3096.604,
    ]
)


@pytest.mark.parametrize("cp_tr_ans", [(4.964700)])
def test_Cp_tr(cp_tr_ans):
    cp_tr_test = cp_tr()
    npt.assert_approx_equal(cp_tr_test, cp_tr_ans, significant=5)


@pytest.mark.parametrize("cp_rot_ans", [(2.978820)])
def test_Cp_rot(cp_rot_ans):
    cp_rot_test = cp_rot()
    npt.assert_approx_equal(cp_rot_test, cp_rot_ans, significant=4)


@pytest.mark.parametrize("freqs, temp, cp_vib_ans", [(ethane_freqs, 298.15, 2.546373)])
def test_Cp_vib(freqs, temp, cp_vib_ans):
    cp_vib_test = cp_vib(freqs, temp)
    npt.assert_approx_equal(cp_vib_test, cp_vib_ans, significant=4)

