#include "md.h"
#include "constants.h"
#include <math.h>

typedef std::vector<std::vector<double>> forcegrid;

void MD::get_accelerations(size_t step)
{
    // Reset our energy-logging data members if this is an output step.
    if (step % d_nstout == 0)
        d_log_LJ_energy = 0;

    // Initialize the forcegrids to 0.
    forcegrid accel_x(d_AtomList.size(), std::vector<double>(d_AtomList.size(), 0));
    forcegrid accel_y(d_AtomList.size(), std::vector<double>(d_AtomList.size(), 0));
    forcegrid accel_z(d_AtomList.size(), std::vector<double>(d_AtomList.size(), 0));

    // Loop over all combinations.
    for (size_t i = 0; i != d_AtomList.size() - 1; ++i)
    {
        for (size_t j = i + 1; j != d_AtomList.size(); ++j)
        {
            // Get the distance in each dimension.
            double r_x = d_AtomList[j].x[0] - d_AtomList[i].x[0];
            double r_y = d_AtomList[j].x[1] - d_AtomList[i].x[1];
            double r_z = d_AtomList[j].x[2] - d_AtomList[i].x[2];

            // Apply minimum image convention.
            if (d_usePBC)
            {
                if (r_x > 0.5 * d_boxsize[0])
                    r_x -= d_boxsize[0];
                else if (r_x <= -0.5 * d_boxsize[0])
                    r_x += d_boxsize[0];

                if (r_y > 0.5 * d_boxsize[1])
                    r_y -= d_boxsize[1];
                else if (r_y <= -0.5 * d_boxsize[1])
                    r_y += d_boxsize[1];

                if (r_z > 0.5 * d_boxsize[2])
                    r_z -= d_boxsize[2];
                else if (r_z <= -0.5 * d_boxsize[2])
                    r_z += d_boxsize[2];
            }

            // Compute the distance.
            double rmag_sq = r_x * r_x + r_y * r_y + r_z * r_z;

            // Compute Lennard-Jones energy and force.
            if (rmag_sq <= d_LJcutoff * d_LJcutoff)
            {
                double rmag = sqrt(rmag_sq);
                
                // Combination rule (arithmetic mean).
                double epsilon = 0.5 * (d_AtomList[i].epsilon + d_AtomList[j].epsilon); // kJ/mol.
                double sigma   = 0.5 * (d_AtomList[i].sigma   + d_AtomList[j].sigma);   // nm.

                // Compute Lennard-Jones force.
                double factor   = (sigma * sigma * sigma * sigma * sigma * sigma) / (rmag_sq * rmag_sq * rmag_sq);
                double attrac   = factor / rmag;
                double repuls   = attrac * factor;
                double LJ_force = 24 * epsilon * (2 * repuls - attrac);

                // Compute Lennard-Jones energy.
                if (step % d_nstout == 0)
                    d_log_LJ_energy += 4 * epsilon * (factor * factor - factor);

                // Decompose force.
                double force_x = LJ_force * (r_x / rmag);
                double force_y = LJ_force * (r_y / rmag);
                double force_z = LJ_force * (r_z / rmag);

                // Fill the force grid (are the indices of the masses correct?).
                accel_x[i][j] =   force_x / d_AtomList[i].mass;
                accel_x[j][i] = - force_x / d_AtomList[j].mass;

                accel_y[i][j] =   force_y / d_AtomList[i].mass;
                accel_y[j][i] = - force_y / d_AtomList[j].mass;

                accel_z[i][j] =   force_z / d_AtomList[i].mass;
                accel_z[j][i] = - force_z / d_AtomList[j].mass;
            }
            else
            {   // Fill the force grid.
                accel_x[i][j] = 0; accel_x[j][i] = 0;
                accel_y[i][j] = 0; accel_y[j][i] = 0;
                accel_z[i][j] = 0; accel_z[j][i] = 0;
            }
        }
    }

    // Reduce forces.
    std::vector<double> reduced_x(d_AtomList.size(), 0);
    std::vector<double> reduced_y(d_AtomList.size(), 0);
    std::vector<double> reduced_z(d_AtomList.size(), 0);
    
    // Parallellizing the reduction of accelerations is thread-safe, 
    // but slower for small systems.
    // #pragma omp parallel for
    for (size_t i = 0; i < d_AtomList.size(); ++i)
    {
        for (size_t j = 0; j < d_AtomList.size(); ++j)
        {
            reduced_x[i] += accel_x[j][i];
            reduced_y[i] += accel_y[j][i];
            reduced_z[i] += accel_z[j][i];
        }
    }

    // Update the accelerations.
    for (size_t i = 0; i != d_AtomList.size(); ++i)
    {
        d_AtomList[i].a_new[0] = reduced_x[i];
        d_AtomList[i].a_new[1] = reduced_y[i];
        d_AtomList[i].a_new[2] = reduced_z[i];
    }
}
