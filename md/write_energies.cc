#include "md.h"
#include "constants.h"

namespace {
// Instantaneous temperature = sum(mv^2) / 2Nk_B
// v is in nm/ps = 1e3 m/s --> v^2 = (nm/ps)^2 = 1e6 (m/s)^2
// m is in amu = 1.66e-27 kg
static double const factor = (1'000'000 * constants::amu) / (2 * constants::kB);
static bool wroteHeader = false;
}

void MD::writeEnergies(size_t step) const
{   // Do loop for temperature and kinetic_energy.
    double mass_vmag_sq = 0;
    for (Atom const &atom : d_AtomList)
    {
        double vmag_sq = atom.v[0] * atom.v[0] + atom.v[1] * atom.v[1] + atom.v[2] * atom.v[2];
        mass_vmag_sq += atom.mass * vmag_sq;
    }
    
    // Compute relevant energies.
    double temperature      = factor * mass_vmag_sq / d_AtomList.size();
    double energy_kinetic   = 0.5 * mass_vmag_sq;
    double work_thermostat  = d_log_thermo_energy;
    double energy_conserved = energy_kinetic - work_thermostat;

    auto file = fopen("energy.log", "a");

    if (not wroteHeader)
    {
        fprintf(file, "step    T      Ekin       Ethermo    Econs\n");
        wroteHeader = true;
    }

    fprintf
    (
        file,
        "%-6zu  %.1f  %.3E  %.3E  %.3E\n", 
        step, 
        temperature,
        energy_kinetic,
        work_thermostat,
        energy_conserved
    );

    fclose(file);
}
