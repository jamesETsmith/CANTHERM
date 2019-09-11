import numpy as np
import pytest

from cantherm.statmech import zpve
from ethane_data import ethane_freqs_all
from pvc_data import pvc_freqs

npt = np.testing


@pytest.mark.parametrize(
    "freqs, temp, zpve_ans",
    [(ethane_freqs_all, 298.15, 0.074378 * 627.509), (pvc_freqs, 298.15, 62.79319)],
)
def test_zpve(freqs, temp, zpve_ans):
    # Answer is from Gaussian and is in Hartrees, zpve is in kcal/mol
    npt.assert_approx_equal(zpve_ans, zpve(freqs, temp, scale=1.0), significant=5)

