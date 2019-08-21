import numpy as np

from cantherm.constants import kb, N_avo, h, c_in_cm, R_kcal, j_to_cal


def H_tr(temp):
    H = 5.0 / 2.0 * R_kcal * temp
    return H


def H_rot(temp):
    H = 3.0 / 2.0 * R_kcal * temp
    return H


def H_vib(freqs, temp, scale=0.99):
    freqs *= scale
    H = 0
    for nu in freqs:
        ei = h * nu * c_in_cm  # hv for this mode in J
        H += ei / (np.exp(ei / (kb * temp)) - 1.0) * (N_avo * j_to_cal / 1e3)
    return H
