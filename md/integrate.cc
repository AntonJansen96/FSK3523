#include "md.h"
#include "constants.h"
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
    double sigma, vmag_sq0, vmag_sq1;

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
        {   
            // Get the velocity of the atom before thermostat.
            vmag_sq0 = atom.v[0] * atom.v[0] + atom.v[1] * atom.v[1] + atom.v[2] * atom.v[2];

            // Standard deviation of Maxwell-Boltzmann distribution.
            sigma = sqrt((constants::kB * d_T) / (atom.mass * constants::amu));
            
            // Generate normal distribution with mean = 0 and sdev = sigma.
            std::normal_distribution<double> distN(0, sigma);

            // Get new velocities from Maxwell-Boltzmann distribution.
            atom.v[0] = factor * distN(d_engine);
            atom.v[1] = factor * distN(d_engine);
            atom.v[2] = factor * distN(d_engine);
            
            // Get Ekin of the atom after thermostat and log energy difference.
            vmag_sq1 = atom.v[0] * atom.v[0] + atom.v[1] * atom.v[1] + atom.v[2] * atom.v[2];
            d_log_thermo_energy += 0.5 * atom.mass * (vmag_sq1 - vmag_sq0);

            #ifdef DEBUG
            printf("integrator: updating velocity of atom %zu using thermostat\n", atom.idx);
            #endif            
        }
        else
        {   // Update velocities like normal.
            atom.v[0] += 0.5 * (atom.a[0] + atom.a_new[0]) * d_dt;
            atom.v[1] += 0.5 * (atom.a[1] + atom.a_new[1]) * d_dt;
            atom.v[2] += 0.5 * (atom.a[2] + atom.a_new[2]) * d_dt;        
            
            #ifdef DEBUG
            printf("integrator: updating velocity of atom %zu like normal:\n", atom.idx);
            auto const a_mag_sq = atom.a[0] * atom.a[0] + atom.a[1] * atom.a[1] + atom.a[2] * atom.a[2];
            printf("integrator: a_mag_sq = %f\n", a_mag_sq);
            #endif
        }

        // Update accelerations.
        std::swap(atom.a, atom.a_new);
        
        #ifdef DEBUG
        printf("integrator: swapping a and a_new\n");
        auto const a_mag_sq = atom.a[0] * atom.a[0] + atom.a[1] * atom.a[1] + atom.a[2] * atom.a[2];
        auto const a_new_mag_sq = atom.a_new[0] * atom.a_new[0] + atom.a_new[1] * atom.a_new[1] + atom.a_new[2] * atom.a_new[2];
        printf("integrator: a_mag_sq = %f, a_new_mag_sq = %f\n", a_mag_sq, a_new_mag_sq);
        #endif
    }
}
