import numpy as np
import pytest

from cantherm.constants import N_avo
from cantherm.core import h_tr, h_rot, h_vib

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


@pytest.mark.parametrize(" T, H_tr_ans", [(298.15, 1.480)])
def test_Htr(T, H_tr_ans):
    H_tr_test = h_tr(T)
    npt.assert_approx_equal(H_tr_test, H_tr_ans, significant=4)


@pytest.mark.parametrize(" T, H_rot_ans", [(298.15, 0.8881)])
def test_Hrot(T, H_rot_ans):
    H_rot_test = h_rot(T)
    npt.assert_approx_equal(H_rot_test, H_rot_ans, significant=4)


@pytest.mark.parametrize(" freqs, T, H_vib_ans", [(ethane_freqs, 298.15, 0.1582)])
def test_Hvib(freqs, T, H_vib_ans):
    H_vib_test = h_vib(freqs, T)
    npt.assert_approx_equal(H_vib_test, H_vib_ans, significant=4)

