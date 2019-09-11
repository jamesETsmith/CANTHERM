import numpy as np
import pytest

from cantherm.constants import N_avo
from cantherm.statmech import s_tr, s_rot, s_vib
from ethane_data import ethane_freqs, ethane_freqs_all, ethane_Iext, ethane_masses
from pvc_data import pvc_masses, pvc_Iext, pvc_freqs

npt = np.testing


@pytest.mark.parametrize(
    "masses, T, s_tr_ans",
    [(ethane_masses, 298.15, 36.109787), (pvc_masses, 298.15, 40.383)],
)
def test_S_tr(masses, T, s_tr_ans):
    s_tr_test = s_tr(masses, T)
    npt.assert_approx_equal(s_tr_test, s_tr_ans, significant=3)


@pytest.mark.parametrize(
    "sigma, I_ext, temp, s_rot_ans",
    [(1, ethane_Iext, 298.15, 19.848), (1, pvc_Iext, 298.15, 28.746)],
)
def test_S_rot(sigma, I_ext, temp, s_rot_ans):
    s_rot_test = s_rot(sigma, I_ext, temp)
    npt.assert_approx_equal(s_rot_test, s_rot_ans, significant=3)


@pytest.mark.parametrize(
    "freqs, temp, scale, s_vib_ans",
    [
        (ethane_freqs.copy(), 298.15, 0.99, 0.6462),
        (ethane_freqs_all, 298.15, 1, 1.996),
        (pvc_freqs.copy(), 298.15, 1.0, 19.427),
    ],
)
def test_S_vib(freqs, temp, scale, s_vib_ans):
    s_vib_test = s_vib(freqs, temp, scale=scale)
    npt.assert_approx_equal(s_vib_test, s_vib_ans, significant=3)

