#!/bin/bash

cd tests

# Ethane Tests (Thermo + IR)
python ethane_test.py

# H-Transf Tests (Thermo + Kin.)
echo
echo
python h-transf_test.py

# PVC Tests (Thermo + IR)
echo
echo
python pvc_test.py

# Programmatic H-Transf Tests (Thermo.)
echo
echo
python programmatic.py
