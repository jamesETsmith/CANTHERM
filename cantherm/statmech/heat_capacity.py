import numpy as np

from scipy.constants import Boltzmann, N_A, h, c, physical_constants, calorie

c_in_cm = c*100
R_cal = physical_constants['molar gas constant'][0]/(calorie)

def cp_tr():
    """Returns the translational contribution to the heat capcity.

    This approximate molecules as an ideal gas particle. For more details see 
    The NIST Reference Database I.D.2 (VII.C.6.) <https://cccbdb.nist.gov/thermo.asp> equation 15
    
    Returns
    -------
    float     
        The translational contribution to the heat capacity at constant pressure in cal/(mol K).
    """
    cp = 5.0 / 2.0 * R_cal
    return cp


def cp_rot():
    """Returns the rotational contribution to the heat capacity. ASSUMES molecules are non-linear.

    This calculation assumes rigid body rotation. Also assumes the molecule is non-linear. For more details see 
    The NIST Reference Database I.D.2 (VII.C.6.) <https://cccbdb.nist.gov/thermo.asp> equation 20
    
    Returns
    -------
    float
        The rotational contribution to the heat capacity at constant pressure in cal/(mol K).
    """
    cp = 5.0 / 2.0 * R_cal  # TODO
    return cp


def cp_vib(freqs, temp, scale=0.99):
    """Calculates the vibrational contribution to the constant pressure heat capacity.

    This calculation assumes vibrations act as harmonic oscillators. For more details see 
    The NIST Reference Database I.D.2 (VII.C.6.) <https://cccbdb.nist.gov/thermo.asp> equation 27
    
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
    freqs = freqs.copy()
    freqs *= scale
    for nu in freqs:
        ei = h * nu * c_in_cm  # hv for this mode in J
        cp += (
            R_cal
            * (ei / (Boltzmann * temp)) ** 2
            * np.exp(ei / (Boltzmann * temp))
            / (1.0 - np.exp(ei / (Boltzmann * temp))) ** 2
        )
    return cp
