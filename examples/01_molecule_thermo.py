from cantherm import get_sample_file_path
from cantherm.chemistry.molecule import CMol

# Boiler plate to point Cantherm toward it's own data directory
geom_opt_file = get_sample_file_path("BF3_geom_opt.log")
uccsd_t_file = get_sample_file_path("BF3_uccsd_t.log")

# Use Cantherm to process the logfiles
cmol_opt = CMol(geom_opt_file)
cmol_energy = CMol(uccsd_t_file)

# Calculate the Gibbs Free Energy for the calculation with frequencies.
# We can calculate the Gibbs correction by G - E, then add it to the accurate energy
sigma = 1  # Check by grep "Rotational Symmetry" ../data/BF3_geom_opt.log
g_opt = cmol_opt.calc_free_energy(298.15, sigma, scale=1.0, units="hartree")
g_correction = g_opt - cmol_opt.energy

print(f"UCCSD(T) Energy            = {cmol_energy.energy:.6f} Ha")
print(f"UCCSD(T) Gibbs Free Energy = {cmol_energy.energy + g_correction:.6f} Ha")

print(cmol_opt.data.entropy)
