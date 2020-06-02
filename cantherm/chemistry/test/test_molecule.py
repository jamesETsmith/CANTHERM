import os
import pytest
import numpy as np
import cclib

import cantherm
from cantherm.chemistry.molecule import CMol

npt = np.testing

# Data
data_files = [
    "bz.log",
    # "dvb_ir.out",
    # "ch4_ccsd_freq_opt.log",
    # "phosphonyl.log",
]  # Devs add new tests here
data_rot_sym = [12, 12, 3]  # Devs add new tests here

data_dir = os.path.join(cantherm.__path__[0], "../data")
data_paths = [os.path.join(data_dir, f) for f in data_files]


def test_load_data():
    for file_path in data_paths:
        exists = os.path.exists(file_path)
        npt.assert_equal(True, exists)


@pytest.mark.parametrize("data_path", data_paths)
def test_cmol_att(data_path):
    # Check mol attributes exists and match cclib data
    cc_data = cclib.io.ccread(data_path)
    cmol = CMol(data_path)

    # Calc prop
    npt.assert_equal(getattr(cmol, "success"), cc_data.metadata["success"])
    npt.assert_equal(getattr(cmol, "basis_set"), cc_data.metadata["basis_set"])
    npt.assert_equal(getattr(cmol, "methods"), cc_data.metadata["methods"])

    # Molecular prop
    npt.assert_equal(True, hasattr(cmol, "energy"))
    npt.assert_equal(True, hasattr(cmol, "masses"))
    npt.assert_equal(True, hasattr(cmol, "mom_inertia"))
    npt.assert_equal(getattr(cmol, "geom"), cc_data.atomcoords[-1])


@pytest.mark.parametrize("data_path, sigma", zip(data_paths, data_rot_sym))
def test_cmol_thermo(data_path, sigma):
    # This test should only include sample outputs with thermo data
    cmol = CMol(data_path)

    # print(cmol.calc_ZPVE(scale=0.99, units="hartree"))
    # print(cmol.calc_entropy(298.15, sigma, scale=1.0))
    print(cmol.calc_free_energy(298.15, sigma, scale=1.0, units="kcal/mol"))

    npt.assert_almost_equal(
        cmol.calc_free_energy(298.15, sigma, scale=1.0, units="hartree"),
        cmol.data.freeenergy,
        decimal=6,
    )

