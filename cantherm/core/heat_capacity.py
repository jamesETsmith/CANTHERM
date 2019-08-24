import numpy as np

from cantherm.constants import kb, N_avo, h, c_in_cm, R_cal


def cp_tr():
    """Returns the translational contribution to the heat capcity.

    This approximate molecules as an ideal gas particle. TODO
    
    Returns
    -------
    float     
        The translational contribution to the heat capacity at constant pressure in cal/(mol K).
    """
    cp = 5.0 / 2.0 * R_cal
    return cp


def cp_rot():
    """Returns the rotational partition function. ASSUMES molecules are non-linear.

    This calculation assumes rigid body rotation. TODO
    
    Returns
    -------
    float
        The rotational contribution to the heat capacity at constant pressure in cal/(mol K).
    """
    cp = 3.0 / 2.0 * R_cal
    return cp


def cp_vib(freqs, temp, scale=0.99):
    """Calculates the vibrational contribution to the constant pressure heat capacity.

    This calculation assumes vibrations act as harmonic oscillators. TODO
    
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
        The vibrational contribution to the heat capacity at constant pressure in cal/(mol K).
    """
    cp = 0
    freqs *= scale
    for nu in freqs:
        ei = h * nu * c_in_cm  # hv for this mode in J
        cp += (
            R_cal
            * (ei / (kb * temp)) ** 2
            * np.exp(ei / (kb * temp))
            / (1.0 - np.exp(ei / (kb * temp))) ** 2
        )
    return cp
