#!/usr/bin/env python

import os
import sys
from Harmonics import *

file = sys.argv[1]


Kcos = []
Ksin = []
harmonic = Harmonics(5, Kcos, Ksin)
harmonic.fitPotential(file)

for i in range(5):
    print(harmonic.Kcos[i], harmonic.Ksin[i])

for i in range(5):
    print(harmonic.Kcos[i], " ", "%3.1f" % (float(i + 1)))
for i in range(5):
    print(harmonic.Ksin[i], " ", "%3.1f" % (float(i + 1)))
