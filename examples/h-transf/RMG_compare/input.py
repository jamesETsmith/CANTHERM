# source activate rmg_env
# cantherm.py input.py

modelChemistry = "CBS-QB3"
frequencyScaleFactor = 0.99
useHinderedRotors = False
useBondCorrections = True

species('GS', 'GS.py')
transitionState('TS', 'TS.py')
thermo('GS', 'NASA')
# thermo('TS', 'NASA')


reaction(
    label = 'Test',
    reactants = ['GS'],
    products = [],
    transitionState = 'TS',
    tunneling='Wigner',
)

statmech('TS')
kinetics(
    label = 'Test',
    Tmin = (298,'K'), Tmax = (1500,'K'), Tcount = 8, # this can be changed to any desired temperature range with any number of temperatures
    Tlist = ([298, 300, 400, 500, 600, 800, 1000, 1500],'K'),
    )
