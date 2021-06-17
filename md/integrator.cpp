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
    double sigma, vmag_sq0;

    for (Atom &atom : d_AtomList)
    {
        // Update positions.
        atom.x[0] += (atom.v[0] + 0.5 * atom.a[0] * d_dt) * d_dt;
        atom.x[1] += (atom.v[1] + 0.5 * atom.a[1] * d_dt) * d_dt;
        atom.x[2] += (atom.v[2] + 0.5 * atom.a[2] * d_dt) * d_dt;

        // Apply periodic boundary condition (PBC).
        if (d_usePBC)
        {
            for (size_t i : {0, 1, 2})
            {
                if (atom.x[i] > d_boxsize[i])
                    atom.x[i] -= d_boxsize[i];
                else if (atom.x[0] < 0.0)
                    atom.x[i] += d_boxsize[i];
            }
        }

        // Update velocities.
        if (d_useThermostat and ((rand() / randMax) < d_tauT * d_dt))
        {   
            // Get the velocity of the atom before thermostat.
            vmag_sq0 = atom.vmag_sq();

            // Standard deviation of Maxwell-Boltzmann distribution.
            sigma = sqrt((constants::kB * d_T) / (atom.mass * constants::amu));
            
            // Generate normal distribution with mean = 0 and sdev = sigma.
            std::normal_distribution<double> distN(0, sigma);

            // Get new velocities from Maxwell-Boltzmann distribution.
            atom.v[0] = factor * distN(d_engine);
            atom.v[1] = factor * distN(d_engine);
            atom.v[2] = factor * distN(d_engine);
            
            // Get Ekin of the atom after thermostat and log energy difference.
            d_log_thermo_energy += 0.5 * atom.mass * (atom.vmag_sq() - vmag_sq0);

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
            printf("integrator: a_mag_sq = %f\n", atom.amag_sq());
            #endif
        }

        // Update accelerations.
        std::swap(atom.a, atom.a_new);
        
        #ifdef DEBUG
        printf("integrator: swapping a and a_new\n");
        printf("integrator: a_mag_sq = %f, a_new_mag_sq = %f\n", atom.amag_sq(), atom.a_new_mag_sq());
        #endif
    }
}
