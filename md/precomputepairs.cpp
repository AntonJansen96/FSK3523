#include "md.h"

void MD::precomputepairs()
{
    // Give the grids the proper size.
    d_pairs_eps.resize(d_AtomList.size(), std::vector<real>(d_AtomList.size()));
    d_pairs_sig.resize(d_AtomList.size(), std::vector<real>(d_AtomList.size()));

    // Loop over all combinations.
    for (size_t i = 0; i != d_AtomList.size() - 1; ++i)
    {
        for (size_t j = i + 1; j != d_AtomList.size(); ++j)
        {
            // For each combination, set the sigma and epsilon.
            d_pairs_eps[i][j] = 0.5 * (d_AtomList[i].epsilon + d_AtomList[j].epsilon);
            d_pairs_sig[i][j] = 0.5 * (d_AtomList[i].sigma   + d_AtomList[j].sigma);
        }
    }
}
