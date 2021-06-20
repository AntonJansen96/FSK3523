#include "md.h"
#include "constants.h"

namespace {
// Instantaneous temperature = sum(mv^2) / 2Nk_B
// v is in nm/ps = 1e3 m/s --> v^2 = (nm/ps)^2 = 1e6 (m/s)^2
// m is in amu = 1.66e-27 kg
static real const factor = (1'000'000 * constants::amu) / (2 * constants::kB);
static real energy_LJ_tail = 0;
static bool wroteHeader   = false;
static bool hasLJtailcorr = false;
}

void MD::writeEnergies(size_t step) const
{   // Compute the Lennard-Jones tail-correction the first time this function is called.
    if (d_useLJ and not hasLJtailcorr)
    {
        energy_LJ_tail = this->tailcorrection();
        hasLJtailcorr  = true;
    }

    // Do loop for temperature and kinetic_energy.
    real mass_vmag_sq = 0;
    for (Atom const &atom : d_AtomList)
        mass_vmag_sq += atom.mass * atom.vmag_sq();
    
    // Compute relevant energies.
    real temperature      = factor * mass_vmag_sq / d_AtomList.size();
    real energy_kinetic   = 0.5 * mass_vmag_sq;
    real work_thermostat  = d_log_thermo_energy;
    real energy_LJ        = d_log_LJ_energy;
    real energy_LJ_total  = energy_LJ + energy_LJ_tail;
    real energy_conserved = (energy_kinetic - work_thermostat) + energy_LJ_total;

    auto file = fopen("energy.log", "a");

    if (not wroteHeader)
    {
        fprintf(file, "step    T      Ekin       Ethermo     E_LJ      E_LJ_corr   E_LJ_tot   Econs\n");
        wroteHeader = true;
    }

    fprintf
    (
        file,
        "%-6zu  %.1f  %.3E  %.3E  %.3E  %.3E  %.3E  %.3E\n", 
        step, 
        temperature,
        energy_kinetic,
        work_thermostat,
        energy_LJ,
        energy_LJ_tail,
        energy_LJ_total,
        energy_conserved
    );

    fclose(file);
}
