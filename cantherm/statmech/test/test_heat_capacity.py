import numpy as np
import pytest

from cantherm.constants import N_A, R_cal
from cantherm.statmech import cp_tr, cp_rot, cp_vib
from ethane_data import ethane_freqs, ethane_freqs_all
from pvc_data import pvc_freqs


npt = np.testing


@pytest.mark.parametrize("cp_tr_ans", [(4.968011), (2.981 + R_cal)])
def test_Cp_tr(cp_tr_ans):
    cp_tr_test = cp_tr()
    npt.assert_approx_equal(cp_tr_test, cp_tr_ans, significant=4)


@pytest.mark.parametrize("cp_rot_ans", [(4.968011), (2.981 + R_cal)])
def test_Cp_rot(cp_rot_ans):
    cp_rot_test = cp_rot()
    npt.assert_approx_equal(cp_rot_test, cp_rot_ans, significant=4)


@pytest.mark.parametrize("freqs, temp, cp_vib_ans", [(ethane_freqs, 298.15, 2.548059)])
def test_Cp_vib(freqs, temp, cp_vib_ans):
    cp_vib_test = cp_vib(freqs, temp)
    npt.assert_approx_equal(cp_vib_test, cp_vib_ans, significant=4)

