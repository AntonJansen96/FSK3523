#include "md.h"
#include "constants.h"
#include <math.h>

void MD::computeforces(size_t step)
{
    // Reset our energy-logging data members if this is an output step.
    if (step % d_nstout == 0)
        d_log_LJ_energy = 0;

    d_Fgrid_x = d_emptyFgrid; // No construction + move, only a copy (= way faster).
    d_Fgrid_y = d_emptyFgrid; // No construction + move, only a copy (= way faster).
    d_Fgrid_z = d_emptyFgrid; // No construction + move, only a copy (= way faster).

    // Loop over all combinations.
    #pragma omp parallel for
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
            double const rmag_sq = r_x * r_x + r_y * r_y + r_z * r_z;

            // Compute Lennard-Jones energy and force.
            if (d_useLJ and rmag_sq <= d_LJcutoff * d_LJcutoff)
            {
                double const rmag = sqrt(rmag_sq);

                // Combination rule (arithmetic mean).
                double const epsilon = d_pairs_eps[i][j];
                double const sigma   = d_pairs_sig[i][j];
                double const sigma6  = sigma * sigma * sigma * sigma * sigma * sigma;

                // Compute Lennard-Jones force.
                double const factor1  = sigma6 / (rmag_sq * rmag_sq * rmag_sq);
                double const attrac   = factor1 / rmag;
                double const repuls   = attrac * factor1;
                double const LJ_force = 24 * epsilon * (2 * repuls - attrac);

                // Compute truncated-shifted Lennard-Jones potential.
                if (step % d_nstout == 0)
                {
                    double factor2 = sigma6 / d_LJcutoff_6;

                    // Lennard-Jones factors for r and rc.
                    double const LJ_factor_r  = factor1 * factor1 - factor1;
                    double const LJ_factor_rc = factor2 * factor2 - factor2;

                    // Add the shifted-trunacted Lennard-Jones energy to total.
                    #pragma omp critical
                    {
                        d_log_LJ_energy += 4 * epsilon * (LJ_factor_r - LJ_factor_rc);
                    }
                }

                // Decompose force.
                double const force_x = LJ_force * (r_x / rmag);
                double const force_y = LJ_force * (r_y / rmag);
                double const force_z = LJ_force * (r_z / rmag);

                // Fill the force grid.
                d_Fgrid_x[i][j] =   force_x;
                d_Fgrid_x[j][i] = - force_x;

                d_Fgrid_y[i][j] =   force_y;
                d_Fgrid_y[j][i] = - force_y;

                d_Fgrid_z[i][j] =   force_z;
                d_Fgrid_z[j][i] = - force_z;
            }
        }
    }

    this->reduceforces();
}
