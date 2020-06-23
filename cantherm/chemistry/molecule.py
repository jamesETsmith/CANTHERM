"""
Bare bones molecule class for Cantherm

Author: James E T Smith <james.smith9113@gmail.com>
Date: 1/8/20
"""

import numpy as np
import cclib
from cclib.parser.utils import convertor
from cclib.method import Nuclear
from cclib.method.nuclear import get_isotopic_masses

from cantherm import statmech
from scipy.constants import N_A, Boltzmann, h, c
from cantherm.statmech.partition_function import q_tr, q_rot, q_vib

la = np.linalg


class CMol:
    def __init__(self, file_path: str):
        self.data = cclib.io.ccread(file_path)
        self.file_path = file_path

        # Calculation properties
        # Check for success of calculation
        metadata = self.data.metadata

        # Make sure the calculation finished successfully
        self.success = metadata["success"]
        if self.success == False:
            raise ValueError(
                "Calculation wasn't successful, we shouldn't process the data"
            )

        # Get basis set
        self.basis_set = metadata["basis_set"]
        # if self._basis_set == "CBSB3":
        #     raise ValueError("CBSQB3 parsing hasn't been implemented yet.")

        # Get methods used (i.e. HF, DFT, CCSD, CCSD-T, etc)
        self.methods = metadata["methods"]
        self.method = self.methods[-1]  # final method used

        # Molecular properties
        if "CC" in self.method:
            self.energy = convertor(self.data.ccenergies[-1], "eV", "hartree")
        else:
            self.energy = convertor(self.data.scfenergies[-1], "eV", "hartree")

        self.geom = self.data.atomcoords[-1]  # Get final geom
        self.masses = get_isotopic_masses(self.data.atomnos)

        nuc = Nuclear(self.data)
        self.mom_inertia = (
            nuc.principal_moments_of_inertia("g_cm_2")[0] * 1e-7
        )  # in kg * m^2
        # ^ Currently need https://github.com/jamesETsmith/cclib/tree/moment_inertia
        # for the above

        # Optional properties
        try:
            self.vibfreqs = self.data.vibfreqs
        except:
            print("No vib data found")

    def calc_ZPVE(self, scale: float = 0.99, units: str = "kcal/mol") -> float:
        zpve = statmech.zpve(self.vibfreqs, scale=scale)
        if units != "kcal/mol":
            zpve = convertor(zpve, "kcal/mol", units)
        return zpve

    def calc_heat_capacity(
        self, temp: float, scale: float = 0.99, units: str = "kcal/mol"
    ) -> float:
        cp = statmech.cp_tr() + statmech.cp_rot()
        cp += statmech.cp_vib(self.vibfreqs, temp, scale=scale)
        if units != "kcal/mol":
            cp = convertor(cp, "kcal/mol", units)
        return cp

    def calc_entropy(
        self, temp: float, sigma: int, scale: float = 0.99, units: str = "kcal/mol"
    ) -> float:
        ent = statmech.s_tr(self.masses, temp)
        ent += statmech.s_rot(sigma, self.mom_inertia, temp)
        ent += statmech.s_vib(self.vibfreqs, temp, scale=scale)
        if units != "kcal/mol":
            ent = convertor(ent, "kcal/mol", units)

        print("TRANS", statmech.s_tr(self.masses, temp))
        print("ROT", statmech.s_rot(sigma, self.mom_inertia, temp))
        print("VIB", statmech.s_vib(self.vibfreqs, temp, scale=scale))
        return ent

    def calc_free_energy(
        self, temp: float, sigma: int, scale: float = 0.99, units: str = "kcal/mol"
    ) -> float:
        free_energy = statmech.g_tr(self.masses, temp)
        free_energy += statmech.g_rot(sigma, self.mom_inertia, temp)
        free_energy += statmech.g_vib(self.vibfreqs, temp, scale=scale)
        free_energy += convertor(self.energy, "hartree", "kcal/mol")
        free_energy += self.calc_ZPVE(scale=scale)
        if units != "kcal/mol":
            free_energy = convertor(free_energy, "kcal/mol", units)
        return free_energy

    def calculate_Q(self, temp: float, sigma: int, I_ext, freqs, scale=0.99):
        """Returns the translational, rotational, and vibrational contribution to the partition functions.

        This approximate molecules as an ideal gas particle. For more details see 
        The NIST Reference Database I.D.2 (VII.C.6.) <https://cccbdb.nist.gov/thermo.asp> equations 12, 13, 17, 25
    
        Returns
        -------
        float
            The translational, rotational, and vibrational contributions to the partition function
        """
        mass = self.masses.sum() / 1e3 / N_A  # Mass in kg / molecule
        Q = q_tr(mass, temp) 
        Q*= q_rot(sigma, I_ext, temp)
        Q*= q_vib(freqs, temp, scale=0.99)

        print(q_tr(mass, temp))
        print(q_rot(sigma, I_ext, temp))
        print(q_vib(freqs,temp,scale=0.99))
        return Q
