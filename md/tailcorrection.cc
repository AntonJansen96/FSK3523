#include "md.h"
#include "constants.h"

double MD::tailcorrection() const
{
    double energy = 0;

    // We loop over all atoms because epsilon and sigma may be different.
    for (size_t i = 0; i != d_AtomList.size(); ++i)
    {
        double epsilon = d_AtomList[i].epsilon;
        double sigma   = d_AtomList[i].sigma;
        
        // Compute the number density rho.
        double volume  = d_boxsize[0] * d_boxsize[1] * d_boxsize[2];
        double rho     = d_AtomList.size() / volume;
        
        // Compute sigma^3 and sigma^9
        double sigma_3 = sigma   * sigma   * sigma;
        double sigma_9 = sigma_3 * sigma_3 * sigma_3;

        // Compute r_c^3 and r_c^9
        double rc_3    = d_LJcutoff * d_LJcutoff * d_LJcutoff;
        double rc_9    = rc_3       * rc_3       * rc_3;

        // Compute pre-factor
        double factor  = (8.0 / 3.0) * constants::pi * rho * epsilon * sigma_3;
        
        // Compute repulsive and attractive terms.
        double repuls  = (1.0 / 3.0) * (sigma_9 / rc_9);
        double attrac  = sigma_3 / rc_3;
        
        energy += factor * (repuls - attrac);
    }
    return energy;
}
