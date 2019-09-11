from cantherm.constants import j_to_cal, c_in_cm, N_avo, h


def zpve(freqs, temp, scale=0.99):
    """Calculated the zero point correction to the energy from vibrations.

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
        The zero point vibrational energy in kcal/mol.
    """
    energy = 0
    freqs = freqs.copy()
    freqs *= scale

    for nu in freqs:
        energy += 0.5 * h * nu * c_in_cm * N_avo * j_to_cal / 1e3  # kcal/mol

    return energy
