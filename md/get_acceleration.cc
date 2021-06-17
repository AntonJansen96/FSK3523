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

void MD::get_accelerations(size_t step)
{
    // Reset our energy-logging data members if this is an output step.
    if (step % d_nstout == 0)
        d_log_LJ_energy = 0;

    // Initialize the forcegrids to 0.
    grid Fgrid_x(d_AtomList.size(), std::vector<double>(d_AtomList.size(), 0));
    grid Fgrid_y(d_AtomList.size(), std::vector<double>(d_AtomList.size(), 0));
    grid Fgrid_z(d_AtomList.size(), std::vector<double>(d_AtomList.size(), 0));

    #ifdef DEBUG
    std::cout << "accelerate: initialized force grids:\n";
    std::cout << "x \n"; printgrid(Fgrid_x);
    std::cout << "y \n"; printgrid(Fgrid_y);
    std::cout << "z \n"; printgrid(Fgrid_z);
    #endif

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
                double epsilon  = d_pairs_epsilon[i][j];
                double sigma    = d_pairs_sigma[i][j];

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
                    #warning "unefficient"
                    double LJ_energy_rc = 4 * epsilon * (pow(sigma / d_LJcutoff, 12) - pow(sigma / d_LJcutoff, 6));

                    // Shifted-trunacted Lennard-Jones energy.
                    d_log_LJ_energy += LJ_energy_r - LJ_energy_rc;

                    #ifdef DEBUG
                    std::cout << "accelerate: LJ-force " << LJ_force << ", LJ-energy " << LJ_energy_r - LJ_energy_rc << '\n';
                    #endif
                }

                // Decompose force.
                double force_x = LJ_force * (r_x / rmag);
                double force_y = LJ_force * (r_y / rmag);
                double force_z = LJ_force * (r_z / rmag);

                // Fill the force grid.
                Fgrid_x[i][j] +=   force_x;
                Fgrid_x[j][i] += - force_x;

                Fgrid_y[i][j] +=   force_y;
                Fgrid_y[j][i] += - force_y;

                Fgrid_z[i][j] +=   force_z;
                Fgrid_z[j][i] += - force_z;

            #ifdef DEBUG
                std::cout << "accelerate: force_x = " << force_x << '\n';
                std::cout << "accelerate: force_y = " << force_y << '\n';
                std::cout << "accelerate: force_z = " << force_z << "\nn";

                std::cout << "accelerate: force_x of " << d_AtomList[i].idx << " on " << d_AtomList[j].idx << " = " << Fgrid_x[i][j] << '\n';
                std::cout << "accelerate: force_x of " << d_AtomList[j].idx << " on " << d_AtomList[i].idx << " = " << Fgrid_x[j][i] << '\n';

                std::cout << "accelerate: force_y of " << d_AtomList[i].idx << " on " << d_AtomList[j].idx << " = " << Fgrid_y[i][j] << '\n';
                std::cout << "accelerate: force_y of " << d_AtomList[j].idx << " on " << d_AtomList[i].idx << " = " << Fgrid_y[j][i] << '\n';

                std::cout << "accelerate: force_z of " << d_AtomList[i].idx << " on " << d_AtomList[j].idx << " = " << Fgrid_z[i][j] << '\n';
                std::cout << "accelerate: force_z of " << d_AtomList[j].idx << " on " << d_AtomList[i].idx << " = " << Fgrid_z[j][i] << '\n';
            }
            else
            {
                std::cout << "accelerate: otuside LJ_cutoff = " << d_LJcutoff << '\n';
            #endif
            }
        }
    }

    // Reduce forces.
    std::vector<double> reduced_x(d_AtomList.size(), 0);
    std::vector<double> reduced_y(d_AtomList.size(), 0);
    std::vector<double> reduced_z(d_AtomList.size(), 0);

    #ifdef DEBUG
    std::cout << "\naccelerate: initialized reduction vectors:\n";
    std::cout << "x \n"; easy::print(reduced_x);
    std::cout << "y \n"; easy::print(reduced_y);
    std::cout << "z \n"; easy::print(reduced_z);

    std::cout << "\naccelerate: the force grids:\n";
    std::cout << "x \n"; printgrid(Fgrid_x);
    std::cout << "y \n"; printgrid(Fgrid_y);
    std::cout << "z \n"; printgrid(Fgrid_z);
    #endif

    // Parallellizing the reduction of accelerations is thread-safe, 
    // but slower for small systems.
    // #pragma omp parallel for
    for (size_t i = 0; i < d_AtomList.size(); ++i)
    {
        for (size_t j = 0; j < d_AtomList.size(); ++j)
        {
            reduced_x[i] += Fgrid_x[j][i];
            reduced_y[i] += Fgrid_y[j][i];
            reduced_z[i] += Fgrid_z[j][i];
        }
    }

    #ifdef DEBUG
    std::cout << "\naccelerate: reduced vectors:\n";
    std::cout << "accelerate: "; std::cout << "x = "; easy::print(reduced_x);
    std::cout << "accelerate: "; std::cout << "y = "; easy::print(reduced_y);
    std::cout << "accelerate: "; std::cout << "z = "; easy::print(reduced_z);
    std::cout << "\naccelerate: update a_new:\n\n";
    #endif

    // Update the accelerations.
    for (size_t i = 0; i != d_AtomList.size(); ++i)
    {
        d_AtomList[i].a_new[0] = reduced_x[i] / d_AtomList[i].mass;
        d_AtomList[i].a_new[1] = reduced_y[i] / d_AtomList[i].mass;
        d_AtomList[i].a_new[2] = reduced_z[i] / d_AtomList[i].mass;

        #ifdef DEBUG
        std::cout << "accelerate: atom " << d_AtomList[i].idx << '\n';
        std::cout << "accelerate: a_new_x   = " << reduced_x[i] / d_AtomList[i].mass << '\n';
        std::cout << "accelerate: a_new_y   = " << reduced_y[i] / d_AtomList[i].mass << '\n';
        std::cout << "accelerate: a_new_z   = " << reduced_z[i] / d_AtomList[i].mass << '\n';
        std::cout << "accelerate: a_new_mag = " << (reduced_x[i] * reduced_x[i] + reduced_y[i] * reduced_y[i] + reduced_z[i] * reduced_z[i]) / d_AtomList[i].mass << "\n\n";
        #endif
    }
}
