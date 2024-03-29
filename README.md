[![Build Status](https://travis-ci.org/jamesETsmith/CANTHERM.svg?branch=v2.0)](https://travis-ci.org/jamesETsmith/CANTHERM)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/6cc3025f25894e988329082059a21c05)](https://www.codacy.com/app/jamesETsmith/CANTHERM?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=jamesETsmith/CANTHERM&amp;utm_campaign=Badge_Grade)
[![codecov](https://codecov.io/gh/jamesETsmith/CANTHERM/branch/v2.0/graph/badge.svg)](https://codecov.io/gh/jamesETsmith/CANTHERM)


# CANTHERM

A **Python3** based software to calculate thermodynamic properties of molecules and rate coefficients of reactions. At the current time Gaussian log files are require to for the hindered rotor scans and frequency calculations and the CanTherm can extract the energy from either Gaussian or MOLPRO files.

---
## Dependencies

- [scipy](https://www.scipy.org/install.html) 
- [numpy](https://numpy.org/install/) 
- [matplotlib](https://matplotlib.org/3.1.1/users/installing.html)
- [cclib>=1.6.3](https://cclib.github.io/how_to_install.html)
- [pytest-cov](https://github.com/pytest-dev/pytest-cov#installation)

---
## Installation

```bash
git clone git@github.com:jamesETsmith/CANTHERM.git
cd CANTHERM
python -m pip install -e .
```

---
## Running the Tests
```bash
python -m pytest -v test/ --cov=cantherm/ 
```
