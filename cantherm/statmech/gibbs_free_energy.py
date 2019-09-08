import numpy as np

from cantherm.statmech import h_tr, h_rot, h_vib, s_tr, s_rot, s_vib
from cantherm.constants import kb, N_avo, h, c_in_cm


def g_tr(masses, temp):
    """Returns the translational contribution to the Gibbs free energy.

    This approximate molecules as an ideal gas particle. 
    
    Parameters
    ----------
    masses : `np.ndarray`
        The masses of the atoms in the molecules. Units should be in g.
    temp : float
        The temperature in K.
    
    Returns
    -------
    float     
        The translational contribution to the Gibbs free energys in kcal/mol.
    """
    g = h_tr(temp) - temp * s_tr(masses, temp) / 1e3
    return g


def g_rot(sigma, I_ext, temp):
    """Returns the rotational contribution to the Gibbs free energy. Not including ZPVE.
    ASSUMES molecules are non-linear.

    This calculation assumes rigid body rotation.
    
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
        The rotational contribution to the Gibbs free energy in kcal/mol.
    """
    g = h_rot(temp) - temp * s_rot(sigma, I_ext, temp) / 1e3
    return g


def g_vib(freqs, temp, scale=0.99):
    """Calculates the vibrational contribution to the Gibbs free energy.

    This calculation assumes vibrations act as harmonic oscillators.
    
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
        The vibrational contribution to the Gibbs free energy in kcal/mol.
    """
    g = h_vib(freqs, temp, scale=scale) - temp * s_vib(freqs, temp, scale=scale) / 1e3
    return g
