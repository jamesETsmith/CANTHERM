atoms = {
    'C': 3,
    'H': 7,
    'O':2,
}

bonds = {}
linear = False
externalSymmetry = 1
spinMultiplicity = 2
opticalIsomers = 1

energy = {
    'CBS-QB3': GaussianLog('../gs.log')
}

geometry = GaussianLog('../gs.log')
frequencies = GaussianLog('../gs.log')
rotors = []
