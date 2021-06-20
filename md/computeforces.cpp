#include "md.h"
#include "constants.h"
#include <math.h>

#ifdef DEBUG
#include <iostream>

namespace {
void printgrid(grid const &input)
{
    for (auto const &vec : input)
        easy::print(vec);
}
} // namespace
#endif

void MD::computeforces(size_t step)
{
    // Reset our energy-logging data members if this is an output step.
    if (step % d_nstout == 0)
        d_log_LJ_energy = 0;

    d_Fgrid_x = d_emptyFgrid; // No construction + move, only a copy (= way faster).
    d_Fgrid_y = d_emptyFgrid; // No construction + move, only a copy (= way faster).
    d_Fgrid_z = d_emptyFgrid; // No construction + move, only a copy (= way faster).

    #ifdef DEBUG
    std::cout << "accelerate: initialized force grids:\n";
    std::cout << "x \n"; printgrid(d_Fgrid_x);
    std::cout << "y \n"; printgrid(d_Fgrid_y);
    std::cout << "z \n"; printgrid(d_Fgrid_z);
    #endif

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
            double rmag_sq = r_x * r_x + r_y * r_y + r_z * r_z;

            #ifdef DEBUG
            std::cout << "\naccelerate: i = " << i << " (atom " << d_AtomList[i].idx << "), j = " << j << " (atom " << d_AtomList[j].idx << ")\n\n";
            std::cout << "accelerate: r_x = " << r_x << '\n';
            std::cout << "accelerate: r_y = " << r_y << '\n';
            std::cout << "accelerate: r_z = " << r_z << '\n';
            std::cout << "accelerate: rmag_sq = " << rmag_sq << ", rmag = " << sqrt(rmag_sq) << '\n';
            #endif

            // Compute Lennard-Jones energy and force.
            if (d_useLJ and rmag_sq <= d_LJcutoff * d_LJcutoff)
            {
                double rmag = sqrt(rmag_sq);
                
                // Combination rule (arithmetic mean).
                double epsilon = d_pairs_eps[i][j];
                double sigma   = d_pairs_sig[i][j];

                #ifdef DEBUG
                std::cout << "accelerate: within LJ_cutoff = " << d_LJcutoff << '\n';
                std::cout << "accelerate: eps = " << epsilon << ", sig = " << sigma << '\n';
                #endif

                // Compute Lennard-Jones force.
                double factor   = (sigma * sigma * sigma * sigma * sigma * sigma) / (rmag_sq * rmag_sq * rmag_sq);
                double attrac   = factor / rmag;
                double repuls   = attrac * factor;
                double LJ_force = 24 * epsilon * (2 * repuls - attrac);

                // Compute truncated-shifted Lennard-Jones potential.
                if (step % d_nstout == 0)
                {
                    // Lennard-Jones potential at r.
                    double LJ_energy_r  = 4 * epsilon * (factor * factor - factor);

                    // Lennard-Jones potential at rc.
                    double LJ_energy_rc = 4 * epsilon * (pow(sigma / d_LJcutoff, 12) - pow(sigma / d_LJcutoff, 6));

                    // Add the shifted-trunacted Lennard-Jones energy to total.
                    #pragma omp critical
                    {
                        d_log_LJ_energy += LJ_energy_r - LJ_energy_rc;
                    }

                    #ifdef DEBUG
                    std::cout << "accelerate: LJ-force " << LJ_force << ", LJ-energy " << LJ_energy_r - LJ_energy_rc << '\n';
                    #endif
                }

                // Decompose force.
                double force_x = LJ_force * (r_x / rmag);
                double force_y = LJ_force * (r_y / rmag);
                double force_z = LJ_force * (r_z / rmag);

                // Fill the force grid.
                d_Fgrid_x[i][j] =   force_x;
                d_Fgrid_x[j][i] = - force_x;

                d_Fgrid_y[i][j] =   force_y;
                d_Fgrid_y[j][i] = - force_y;

                d_Fgrid_z[i][j] =   force_z;
                d_Fgrid_z[j][i] = - force_z;

            #ifdef DEBUG
                std::cout << "accelerate: force_x = " << force_x << '\n';
                std::cout << "accelerate: force_y = " << force_y << '\n';
                std::cout << "accelerate: force_z = " << force_z << "\nn";

                std::cout << "accelerate: force_x of " << d_AtomList[i].idx << " on " << d_AtomList[j].idx << " = " << d_Fgrid_x[i][j] << '\n';
                std::cout << "accelerate: force_x of " << d_AtomList[j].idx << " on " << d_AtomList[i].idx << " = " << d_Fgrid_x[j][i] << '\n';

                std::cout << "accelerate: force_y of " << d_AtomList[i].idx << " on " << d_AtomList[j].idx << " = " << d_Fgrid_y[i][j] << '\n';
                std::cout << "accelerate: force_y of " << d_AtomList[j].idx << " on " << d_AtomList[i].idx << " = " << d_Fgrid_y[j][i] << '\n';

                std::cout << "accelerate: force_z of " << d_AtomList[i].idx << " on " << d_AtomList[j].idx << " = " << d_Fgrid_z[i][j] << '\n';
                std::cout << "accelerate: force_z of " << d_AtomList[j].idx << " on " << d_AtomList[i].idx << " = " << d_Fgrid_z[j][i] << '\n';
            #endif
            }
        }
    }

    this->reduceforces();
}
