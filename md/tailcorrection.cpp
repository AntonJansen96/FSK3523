#include "md.h"
#include "constants.h"

real MD::tailcorrection() const
{
    real energy = 0;

    // Compute the number density rho outside of loop.
    real const volume = d_boxsize[0] * d_boxsize[1] * d_boxsize[2];
    real const rho    = d_AtomList.size() / volume;

    // We loop over all atoms because epsilon and sigma may be different.
    for (Atom const &atom : d_AtomList)
    {   // Compute sigma^3 and sigma^9
        real sigma_3 = atom.sigma * atom.sigma * atom.sigma;
        real sigma_9 = sigma_3    * sigma_3    * sigma_3;

        // Compute r_c^3 and r_c^9
        real rc_3    = d_LJcutoff * d_LJcutoff * d_LJcutoff;
        real rc_9    = rc_3       * rc_3       * rc_3;

        // Compute pre-factor
        real factor  = (8.0 / 3.0) * constants::pi * rho * atom.epsilon * sigma_3;
        
        // Compute repulsive and attractive terms.
        real repuls  = (1.0 / 3.0) * (sigma_9 / rc_9);
        real attrac  = sigma_3 / rc_3;
        
        energy += factor * (repuls - attrac);
    }
    return energy;
}
