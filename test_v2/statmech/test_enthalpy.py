import numpy as np
import pytest

from cantherm.constants import N_avo
from cantherm.statmech import h_tr, h_rot, h_vib, zpve
from ethane_data import ethane_freqs
from pvc_data import pvc_freqs

npt = np.testing


@pytest.mark.parametrize(" T, H_tr_ans", [(298.15, 1.481212)])
def test_Htr(T, H_tr_ans):
    H_tr_test = h_tr(T)
    npt.assert_approx_equal(H_tr_test, H_tr_ans, significant=4)


@pytest.mark.parametrize(" T, H_rot_ans", [(298.15, 0.888727)])
def test_Hrot(T, H_rot_ans):
    H_rot_test = h_rot(T)
    npt.assert_approx_equal(H_rot_test, H_rot_ans, significant=4)


@pytest.mark.parametrize(" freqs, T, H_vib_ans", [(ethane_freqs, 298.15, 0.1582)])
def test_Hvib(freqs, T, H_vib_ans):
    H_vib_test = h_vib(freqs, T)
    npt.assert_approx_equal(H_vib_test, H_vib_ans, significant=4)


@pytest.mark.parametrize(
    "freqs, scale, temp, H_tot_ans", [(pvc_freqs, 1.0, 298.15, 0.108739)]
)
def test_Htot(freqs, scale, temp, H_tot_ans):
    H_tot = (
        h_tr(temp)
        + h_rot(temp)
        + h_vib(freqs, temp, scale=scale)
        + zpve(freqs, temp, scale=scale)
    ) / 627.509  # convert to hartree
    npt.assert_approx_equal(H_tot_ans, H_tot, significant=5)
