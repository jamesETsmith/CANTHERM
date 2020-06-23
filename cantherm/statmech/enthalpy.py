import numpy as np
from scipy.constants import Boltzmann, N_A, h, c, calorie, physical_constants

c_in_cm = c*100
R_kcal = physical_constants['molar gas constant'][0]/(calorie*1e3)

def h_tr(temp):
    """Calculates the translational enthalpic contribution.

    Assumes the molecule is an ideal gas. For more details see
    The NIST Reference Database I.D.2 (VII.C.6.) <https://cccbdb.nist.gov/thermo.asp> equation 16
    
    Parameters
    ----------
    temp : float
        The temperature in K.
    
    Returns
    -------
    float
        The translational enthalpy contribution in kcal/mol.
    """
    H = 5.0 / 2.0 * R_kcal * temp
    return H


def h_rot(temp):
    """Calculates the rotational enthalpic contribution.

    Assumes the molecule is an ideal gas AND non-linear. For more details see
    The NIST Reference Database I.D.2 (VII.C.6.) <https://cccbdb.nist.gov/thermo.asp> equation 21
    
    Parameters
    ----------
    temp : float
        The temperature in K.
    
    Returns
    -------
    float
        The rotational enthalpy contribution in kcal/mol.
    """
    H = 3.0 / 2.0 * R_kcal * temp
    return H


def h_vib(freqs, temp, scale=0.99):
    """Calculates the vibrational enthalpic contribution.

    For more details see
    The NIST Reference Database I.D.2 (VII.C.6.) <https://cccbdb.nist.gov/thermo.asp> equation 28
    
    Parameters
    ----------
    freqs : iterable (list or `np.ndarray`)
        A list of the vibrational frequencies in cm^-1.
    temp : float
        The temperature in K.
    scale : float, optional
        Scale of the frequencies, by default 0.99
    
    Returns
    -------
    float 
        The vibrational enthalpy contribution in kcal/mol
    """
    freqs = freqs.copy()
    freqs *= scale
    H = 0
    for nu in freqs:
        ei = h * nu * c_in_cm  # hv for this mode in J
        H += ei / (np.exp(ei / (Boltzmann * temp)) - 1.0) * (N_A / (calorie*1e3))
    return H
