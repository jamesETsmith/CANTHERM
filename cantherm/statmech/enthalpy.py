import numpy as np

from cantherm.constants import kb, N_avo, h, c_in_cm, R_kcal, j_to_cal


def h_tr(temp):
    """Calculates the translational enthalpic contribution.

    Assumes the molecule is an ideal gas.
    
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

    Assumes the molecule is an ideal gas AND non-linear.
    
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
        H += ei / (np.exp(ei / (kb * temp)) - 1.0) * (N_avo * j_to_cal / 1e3)
    return H
