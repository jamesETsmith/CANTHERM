import numpy as np
import pytest

from cantherm.constants import N_avo
from cantherm.statmech import s_tr, s_rot, s_vib

npt = np.testing

ethane_masses = np.array(
    [12.00000, 1.00783, 1.00783, 1.00783, 12.00000, 1.00783, 1.00783, 1.00783]
)
ethane_Iext = np.array([6.272, 25.377, 25.377]) / (1e23 * N_avo)
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

ethane_freqs_all = np.array(
    [
        307.7063,
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


@pytest.mark.parametrize("masses, T, s_tr_ans", [(ethane_masses, 298.15, 36.110)])
def test_S_tr(masses, T, s_tr_ans):
    s_tr_test = s_tr(masses, T)
    npt.assert_approx_equal(s_tr_test, s_tr_ans, significant=5)


@pytest.mark.parametrize(
    "sigma, I_ext, temp, s_rot_ans", [(1, ethane_Iext, 298.15, 19.848)]
)
def test_S_rot(sigma, I_ext, temp, s_rot_ans):
    s_rot_test = s_rot(sigma, I_ext, temp)
    npt.assert_approx_equal(s_rot_test, s_rot_ans, significant=3)


@pytest.mark.parametrize(
    "freqs, temp, scale, s_vib_ans",
    [(ethane_freqs, 298.15, 0.99, 0.6462), (ethane_freqs_all, 298.15, 1, 1.996)],
)
def test_S_vib(freqs, temp, scale, s_vib_ans):
    s_vib_test = s_vib(freqs, temp, scale=scale)
    npt.assert_approx_equal(s_vib_test, s_vib_ans, significant=3)

