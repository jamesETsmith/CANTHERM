import numpy as np
import pytest

from cantherm.statmech import g_tr, g_rot, g_vib
from ethane_data import ethane_freqs, ethane_Iext, ethane_masses, ethane_freqs_all

npt = np.testing


@pytest.mark.parametrize(
    "masses, temp, g_tr_ans",
    [(ethane_masses, 298.15, (1.480225 - 298.15 * 36.109787 / 1e3))],
)
def test_Gtr(masses, temp, g_tr_ans):
    npt.assert_approx_equal(g_tr(masses, temp), g_tr_ans, significant=5)


@pytest.mark.parametrize(
    "sigma, I_ext, temp, g_rot_ans",
    [(6, ethane_Iext, 298.15, 0.888135 - 298.15 * 16.276204 / 1e3)],
)
def test_Grot(sigma, I_ext, temp, g_rot_ans):
    npt.assert_approx_equal(g_rot(sigma, I_ext, temp), g_rot_ans, significant=5)


@pytest.mark.parametrize(
    "freqs, temp, g_vib_ans",
    [(ethane_freqs, 298.15, 0.158189 - 298.15 * 0.646233 / 1e3)],
)
def test_Gvib(freqs, temp, g_vib_ans):
    npt.assert_approx_equal(g_vib(freqs, temp, scale=0.99), g_vib_ans, significant=4)


@pytest.mark.parametrize(
    "masses, sigma, I_ext, freqs, temp, g_tot_ans",
    [
        (
            ethane_masses,
            1,
            ethane_Iext,
            ethane_freqs_all,
            298.15,
            (0.051261 - 0.074378) * 627.509,
        )
    ],
)
def test_Gtot(masses, sigma, I_ext, freqs, temp, g_tot_ans):
    g = g_tr(masses, temp) + g_rot(sigma, I_ext, temp) + g_vib(freqs, temp, scale=1.0)
    npt.assert_approx_equal(g, g_tot_ans, significant=4)

