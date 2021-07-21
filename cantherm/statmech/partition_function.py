import numpy as np

from scipy.constants import Boltzmann, N_A, h, c

c_in_cm = c * 100


def q_tr(masses, temp):
    """Returns the translational contribution to the partition functions.

    This approximate molecules as an ideal gas particle. For more details see 
    The NIST Reference Database I.D.2 (VII.C.6.) <https://cccbdb.nist.gov/thermo.asp> equation 12 and 13
    or
    "Molecular Driving Forces" by Dill equation 11.18 and 11.19.
    
    Parameters
    ----------
    masses : `np.ndarray`
        The masses of the atoms in the molecules. Units should be in g.
    temp : float
        The temperature in K.
    
    Returns
    -------
    float     
        The translational contribution to the partition function.
    """
    # Translational Contrib.
    # TODO Assuming Unimolecular for now so it's technically per volume
    mass = masses.sum() / 1e3 / N_A  # Mass in kg / molecule
    q = ((2 * np.pi * mass) / h ** 2) ** 1.5 / 101325 * (Boltzmann * temp) ** 2.5
    return q


def q_rot(sigma, I_ext, temp):
    """Returns the rotational partition function. ASSUMES molecules are non-linear.

    This calculation assumes rigid body rotation. For more details see 
    The NIST Reference Database I.D.2 (VII.C.6.) <https://cccbdb.nist.gov/thermo.asp> equation 17
    or
    "Molecular Driving Forces" by Dill equation 11.31.
    
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
        The rotational partition function.
    """
    q = 1
    if np.min(I_ext) == 0:  # linear
        q *= (8 * np.pi ** 2 * np.max(I_ext) * Boltzmann * temp) / (sigma * h ** 2)
    else:  # non-linear
        q *= np.power(np.pi * I_ext[0] * I_ext[1] * I_ext[2], 0.5) / sigma
        q *= np.power(8.0 * np.pi ** 2 * Boltzmann * temp / h ** 2, 1.5)
    return q


def q_vib(freqs, temp, scale=0.99):
    """Calculates the vibrational contribution to the partition function.

    This calculation assumes vibrations act as harmonic oscillators. For more details see 
    The NIST Reference Database I.D.2 (VII.C.6.) <https://cccbdb.nist.gov/thermo.asp> equation 25
    or
    "Molecular Driving Forces" by Dill equation 11.26.
    
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
        The vibrational contribution to the partition function.
    """
    q = 1.0
    freqs = freqs.copy()
    freqs *= scale
    for nu in freqs:
        ei = h * nu * c_in_cm  # hv for this mode in J
        q *= 1.0 / (1.0 - np.exp(-ei / (Boltzmann * temp)))
    return q
