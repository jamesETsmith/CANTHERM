import numpy as np
import pytest

from cantherm.constants import N_avo
from cantherm.statmech import q_tr, q_rot, q_vib

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


@pytest.mark.parametrize("masses, T, q_tr_ans", [(ethane_masses, 298.15, 0.647373e7)])
def test_Qtr(masses, T, q_tr_ans):
    q_tr_test = q_tr(masses, T)
    npt.assert_approx_equal(q_tr_test, q_tr_ans, significant=5)


@pytest.mark.parametrize(
    "sigma, I_ext, temp, q_rot_ans", [(1, ethane_Iext, 298.15, 0.485478e4)]
)
def test_Qrot(sigma, I_ext, temp, q_rot_ans):
    q_rot_test = q_rot(sigma, I_ext, temp)
    npt.assert_approx_equal(q_rot_test, q_rot_ans, significant=4)


@pytest.mark.parametrize(
    "freqs, temp, scale, q_vib_ans",
    [(ethane_freqs, 298.15, 0.99, 1.060162), (ethane_freqs_all, 298.15, 1.0, 1.36705)],
)
def test_Q_vib(freqs, temp, scale, q_vib_ans):
    q_vib_test = q_vib(freqs, temp, scale=scale)
    npt.assert_approx_equal(q_vib_test, q_vib_ans, significant=5)
