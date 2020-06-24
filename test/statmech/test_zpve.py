import numpy as np
import pytest

from cantherm.statmech import zpve
from ethane_data import ethane_freqs_all
from pvc_data import pvc_freqs

npt = np.testing


@pytest.mark.parametrize(
    "freqs, zpve_ans", [(ethane_freqs_all, 0.074378 * 627.509), (pvc_freqs, 62.79319)]
)
def test_zpve(freqs, zpve_ans):
    # Answer is from Gaussian and is in Hartrees, zpve is in kcal/mol
    npt.assert_approx_equal(zpve_ans, zpve(freqs, scale=1.0), significant=5)

