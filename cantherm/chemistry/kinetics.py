import math, re
from scipy.constants import h, c, Boltzmann, calorie, physical_constants, N_A
c_in_cm = c*100
R_kcal = physical_constants['molar gas constant'][0]/(calorie*1e3)
ha_to_kcal = physical_constants['Hartree energy'][0]*N_A/(calorie*1e3)

class Reaction:
    """  
    Reaction class for CANTHERM

    Author: Virginia Johnson <virginia.johnson@colorado.edu>
    Date: 1/26/2020
    """

    def __init__(self,
                 reactants,
                 ts,
                 temp,
                 reac_type='Unimol',
                 products=[],
                 tunneling=None,
                 scale=0.99):
        self.reactants = reactants
        self.ts = ts
        self.temp = temp     #YES
        self.reac_type = reac_type     #YES
        self.products = products
        self.tunneling = tunneling
        self.scale = scale
        self.rates = None     #YES
        self.tunneling_coeff = None     #YES
        self.q_ratio = None     #YES
        return

    ############################################################################

    def calc_TST_rates(self):
        """Calculates the transition state theory rate constants for the
            reaction at the given temperatures.
        """
        print('hello0')
        self.rates = [0] * len(self.temp)
        self.tunneling_coeff = [0] * len(self.temp)
        self.q_ratio = [0] * len(self.temp)

        print('hello1')

        for i in range(len(self.temp)):
            t = self.temp[i]
            if self.reac_type == 'Unimol':
                Q_react = self.reactants.calculate_Q(t)
                Q_TS = self.ts.calculate_Q(t)

            self.q_ratio[i] = (Boltzmann * t / h) * (Q_TS / Q_react)
            self.rates[i] = (Boltzmann * t / h) * (Q_TS / Q_react)
            self.rates[i] *= math.exp(-(self.ts.Energy - \
                                        self.reactants.Energy) * ha_to_kcal \
                                        / R_kcal / t)

            if self.tunneling == "Wigner":
                kappa = wigner_correction(t, self.ts.imagFreq, self.scale)
                self.rates[i] *= kappa
                self.tunneling_coeff[i] = kappa

            #print(self.q_ratio[i])
            print('hello2')


################################################################################


def wigner_correction(t, freq, scale):
    # see doi:10.1103/PhysRev.40.749 and doi:10.1039/TF9595500001
    return (1.0 + 1.0 / 24.0 * (h * abs(freq) * scale * c_in_cm / (t * Boltzmann))**2)

