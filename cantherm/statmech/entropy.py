import numpy as np

from scipy.constants import Boltzmann, N_A, h, c, calorie, physical_constants
from cantherm.statmech import q_tr, q_rot, q_vib, h_rot

c_in_cm = c*100
R_cal = physical_constants['molar gas constant'][0]/(calorie)
R_kcal = R_cal/1e3


def s_tr(masses, temp):
    """Calculates the translational entropic contribution.

    Assumes the molecule is an ideal gas. Sackur-Tetrode equation for more details, or see 
    The NIST Reference Database I.D.2 (VII.C.6.) <https://cccbdb.nist.gov/thermo.asp> equation 14
    
    Parameters
    ----------
    masses : `np.ndarray`
        The masses of the atoms in the molecules. Units should be in g.
    temp : float
        The temperature in K.
    
    Returns
    -------
    float
        The translational entropy contribution in cal/(mol K).
    """
    s = R_cal * (np.log(q_tr(masses, temp)) + 5.0 / 2.0)
    return s


def s_rot(sigma, I_ext, temp):
    """Calculates the rotational entropic contribution.

    Assumes the molecule is an ideal gas AND non-linear. For more details see
    The NIST Reference Database I.D.2 (VII.C.6.) <https://cccbdb.nist.gov/thermo.asp> equation 19
    
    Parameters
    ----------
    sigma : int
        The rotational symmetry factor of the molecule.
    I_ext : `iterable` (list or `np.ndarray`)
        The moments of interia for the molecule. Units should be amu * ang.^2.
    temp : float
        The temperature in K.
    
    Returns
    -------
    float
        The rotational entropic contribution in cal/(mol K).
    """
    s = h_rot(temp) * 1e3 / temp + R_cal * np.log(q_rot(sigma, I_ext, temp))
    return s


def s_vib(freqs, temp, scale=0.99):
    """Calculates the vibrational entropic contribution.

   For more details see
   The NIST Reference Database I.D.2 (VII.C.6.) <https://cccbdb.nist.gov/thermo.asp> equation 26
    
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
        The vibrational entropic contribution in cal/(mol K)
    """
    freqs = freqs.copy()
    freqs *= scale
    s = 0
    for nu in freqs:
        ei = h * nu * c_in_cm  # hv for this mode in J
        s += R_cal * (
            (ei / (Boltzmann * temp)) / (np.exp(ei / (Boltzmann * temp)) - 1.0)
            - np.log(1.0 - np.exp(-ei / (Boltzmann * temp)))
        )

    return s
