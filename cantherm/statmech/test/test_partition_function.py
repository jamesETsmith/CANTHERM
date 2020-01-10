import numpy as np
import pytest

from cantherm.constants import N_avo
from cantherm.statmech import q_tr, q_rot, q_vib
from ethane_data import ethane_masses, ethane_Iext, ethane_freqs, ethane_freqs_all
from pvc_data import pvc_masses, pvc_Iext, pvc_freqs

npt = np.testing


@pytest.mark.parametrize(
    "masses, T, q_tr_ans",
    [(ethane_masses, 298.15, 0.647373e7), (pvc_masses, 298.15, 0.549262e08)],
)
def test_Qtr(masses, T, q_tr_ans):
    q_tr_test = q_tr(masses, T)
    npt.assert_approx_equal(q_tr_test, q_tr_ans, significant=5)


@pytest.mark.parametrize(
    "sigma, I_ext, temp, q_rot_ans",
    [(1, ethane_Iext, 298.15, 0.485478e4), (1, pvc_Iext, 298.15, 0.427524e06)],
)
def test_Qrot(sigma, I_ext, temp, q_rot_ans):
    q_rot_test = q_rot(sigma, I_ext, temp)
    npt.assert_approx_equal(q_rot_test, q_rot_ans, significant=4)


@pytest.mark.parametrize(
    "freqs, temp, scale, q_vib_ans",
    [
        (ethane_freqs.copy(), 298.15, 0.99, 1.060162),
        (ethane_freqs_all, 298.15, 1.0, 1.36705),
        (pvc_freqs, 298.15, 1.0, 0.986612e02),
    ],
)
def test_Q_vib(freqs, temp, scale, q_vib_ans):
    q_vib_test = q_vib(freqs, temp, scale=scale)
    npt.assert_approx_equal(q_vib_test, q_vib_ans, significant=5)
