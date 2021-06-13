#include "md.h"

namespace {
// Instantaneous temperature = sum(mv^2) / 2Nk_B
// v is in nm/ps = 1e3 m/s --> v^2 = (nm/ps)^2 = 1e6 (m/s)^2
// m is in amu = 1.66e-27 kg
static double const factor = (1'000'000 * constants::amu) / (2 * constants::kB);
}

void MD::writeEnergies(size_t step)
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

    fprintf(file, "%-6d  %.1f  %.3E\n", step, temp, Ekin);

    fclose(file);
}
