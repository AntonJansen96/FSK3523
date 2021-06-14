#include "md.h"

namespace {
// Instantaneous temperature = sum(mv^2) / 2Nk_B
// v is in nm/ps = 1e3 m/s --> v^2 = (nm/ps)^2 = 1e6 (m/s)^2
// m is in amu = 1.66e-27 kg
static double const factor = (1'000'000 * constants::amu) / (2 * constants::kB);
static bool wroteHeader = false;
}

void MD::writeEnergies(size_t step) const
{
    double vmag_sq, mass_vmag_sq = 0;

    for (Atom const &atom : d_AtomList)
    {
        vmag_sq = atom.v[0] * atom.v[0] + atom.v[1] * atom.v[1] + atom.v[2] * atom.v[2];
        
        mass_vmag_sq += atom.mass * vmag_sq;
    }
    
    double temp = factor * mass_vmag_sq / d_AtomList.size();
    double Ekin = 0.5 * mass_vmag_sq;
    
    auto file = fopen("energy.log", "a");

    if (not wroteHeader)
    {
        fprintf(file, "step    T      Ekin       LJ-trunc    LJ-tail     LJ-total\n");
        wroteHeader = true;
    }

    fprintf
    (
        file,
        "%-6zu  %.1f  %.3E  %.3E  %.3E  %.3E\n", 
        step, temp, Ekin, 
        d_log_LJ_trun_energy, 
        d_log_LJ_tail_energy, 
        d_log_LJ_trun_energy + d_log_LJ_tail_energy
    );

    fclose(file);
}
