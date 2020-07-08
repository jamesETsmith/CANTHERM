import os
import pytest
import numpy as np
import cclib

import cantherm
from cantherm.chemistry.molecule import CMol
from cantherm import get_sample_file_path

npt = np.testing

# Data
data_files = [
    "bz.log",  # Rot. symm. 12
    "BF3_geom_opt.log",  # Rot symm. 1
    # "dvb_ir.out",
    "ch4_ccsd_freq_opt.log",  # Rot. Symm 12
    "phosphonyl.log",  # Rot. Symm 3
    "BF3_freq_orca.out",  # Rot. Symm 3
    "OH_opt_orca.out",  # Rot symm. 2
    "oh_freq.log",  # Rot. symm. 1
    "CuCl.log",  # Rot. symm. 1
]  # Devs add new tests here
data_rot_sym = [12, 1, 12, 3, 3, 1, 1, 1]  # Devs add new tests here

if len(data_files) != len(data_rot_sym):
    raise ValueError(
        "data_files and data_rot_sym are different length!! Edit test_molecule.py"
    )

data_paths = [get_sample_file_path(f) for f in data_files]


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
    npt.assert_equal(getattr(cmol, "basis_set"), cc_data.metadata.get("basis_set", ""))
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

    ZPVE = cmol.calc_ZPVE(scale=1.0, units="hartree")
    entropy = cmol.calc_entropy(298.15, sigma, scale=1.0, units="hartree")
    # The ORCA entropy is weird and doesn't have ZPVE
    if cmol.data.metadata["package"] != "ORCA":
        npt.assert_almost_equal(entropy, cmol.data.entropy, decimal=5)
        npt.assert_almost_equal(ZPVE, cmol.data.zpve, decimal=6)

    err = (
        cmol.calc_free_energy(298.15, sigma, scale=1.0, units="hartree")
        - cmol.data.freeenergy
    )

    # This is loose because there is disagreement with G16 on the mHa level
    if abs(err) > 1.9e-3:
        print(abs(err))
        raise AssertionError("Cantherm's free energy doesn't match the known value")
