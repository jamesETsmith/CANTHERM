# Outline for v2.0

## Stat. Mech. Library 
- [  ] Calculation of partition functions
    - [X] Translational
    - [X] Vibrational
    - [X] Rotational
    - [ ] Hindered Vibrations

- [  ] Calculation of enthalpy
    - [X] Translational
    - [X] Vibrational
    - [X] Rotational
    - [ ] Hindered Vibrations

- [  ] Calculation of entropy
    - [X] Translational
    - [X] Vibrational
    - [X] Rotational
    - [ ] Hindered Vibrations

- [  ] Calculation of gibbs free energy
    - [X] Translational
    - [X] Vibrational
    - [X] Rotational
    - [ ] Hindered Vibrations

- [  ] Calculation of heat capacity at constant temperature
    - [X] Translational
    - [X] Vibrational
    - [X] Rotational
    - [ ] Hindered Vibrations

- [X] Change name to statmech
- [X] Make all functions that use freqs copy
  
## Kinetics
- [ ] Reaction
- [ ] Reaction system
- [ ] 

### Utils
- [ ] Get rid of defined constants, try to convert as many as possible to those from scipy

## Interfaces
- [ ] Use `cclib` to read (most) QC software formats for those that don't have python interface, PySCF for example. The packages with a Python interface should package the necessary data in numpy arrays.
- [ ] Compare PySF to Gaussian

