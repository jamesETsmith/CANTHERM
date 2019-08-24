import numpy as np
import pytest

from cantherm.core import zpve

npt = np.testing

ethane_freqs = np.array(
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


@pytest.mark.parametrize("freqs, temp, zpve_ans", [(ethane_freqs, 298.15, 0.074378)])
def test_zpve(freqs, temp, zpve_ans):
    # Answer is from Gaussian and is in Hartrees, zpve is in kcal/mol
    npt.assert_approx_equal(
        zpve_ans, zpve(freqs, temp, scale=1.0) / 627.509, significant=5
    )

