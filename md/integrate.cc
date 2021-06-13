#include "md.h"
#include <random>
#include <math.h>

namespace {
// We go from m/s to nm/ps so we have a factor 0.001. Furthermore, we have 
// a sqrt(2) from 0.5mv^2 and 1/sqrt(3) from vector calculus.
static double const factor  = 0.001 * sqrt(2.0 / 3.0);
// Do not do this static cast in the loop.
static double const randMax = static_cast<double>(RAND_MAX);
} // Namespace.

void MD::integrate()
{
    for (Atom &atom : d_AtomList)
    {
        // Update positions.
        atom.x[0] += (atom.v[0] + 0.5 * atom.a[0] * d_dt) * d_dt;
        atom.x[1] += (atom.v[1] + 0.5 * atom.a[1] * d_dt) * d_dt;
        atom.x[2] += (atom.v[2] + 0.5 * atom.a[2] * d_dt) * d_dt;

        // Apply periodic boundary condition.
        if (d_usePBC)
        {
            if (atom.x[0] >  d_boxsize[0])
                atom.x[0] -= d_boxsize[0];
            else if (atom.x[0] < 0)
                atom.x[0] += d_boxsize[0];

            if (atom.x[1] >  d_boxsize[1])
                atom.x[1] -= d_boxsize[1];
            else if (atom.x[1] < 0)
                atom.x[1] += d_boxsize[1];

            if (atom.x[2] >  d_boxsize[2])
                atom.x[2] -= d_boxsize[2];
            else if (atom.x[2] < 0)
                atom.x[2] += d_boxsize[2];
        }

        // Update velocities.
        if (d_useThermostat and ((rand() / randMax) < d_tauT * d_dt))
        {   // If we use the Andersen thermostat...
            // Standard deviation of Maxwell-Boltzmann distribution.
            double sigma = sqrt((constants::kB * d_T) / (atom.mass * constants::amu));
            
            // Generate normal distribution with mean = 0 and sdev = sigma.
            std::normal_distribution<double> distN(0, sigma);

            // Generate velocities from Maxwell-Boltzmann distribution.
            atom.v[0] = factor * distN(d_engine);
            atom.v[1] = factor * distN(d_engine);
            atom.v[2] = factor * distN(d_engine);
        }
        else
        {   // Update velocities like normal.
            atom.v[0] += 0.5 * (atom.a[0] + atom.a_new[0]) * d_dt;
            atom.v[1] += 0.5 * (atom.a[1] + atom.a_new[1]) * d_dt;
            atom.v[2] += 0.5 * (atom.a[2] + atom.a_new[2]) * d_dt;
        }

        // Update accelerations.
        std::swap(atom.a, atom.a_new);
    }
}
