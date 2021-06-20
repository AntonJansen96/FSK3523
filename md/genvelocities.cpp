#include "md.h"
#include "constants.h"
#include <math.h>

namespace {
// We go from m/s to nm/ps so we have a factor 0.001. Furthermore, we have 
// a sqrt(2) from 0.5mv^2 and 1/sqrt(3) from vector calculus.
static real const factor = 0.001 * sqrt(2.0 / 3.0);
} // Namespace.

void MD::generate_velocities()
{
    // Loop through the atoms.
    for (Atom &atom : d_AtomList)
    {
        // Standard deviation of Maxwell-Boltzmann distribution.
        real sigma = sqrt((constants::kB * d_T) / (atom.mass * constants::amu));
        
        // Generate normal distribution with mean = 0 and sdev = sigma.
        std::normal_distribution<real> distN(0, sigma);

        // Generate velocities.
        atom.v[0] = factor * distN(d_engine);
        atom.v[1] = factor * distN(d_engine);
        atom.v[2] = factor * distN(d_engine);
    }
}
