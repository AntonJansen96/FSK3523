#include "md.h"
#include <math.h>

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
            // e_ij = sqrt(e_i * e_j)
            d_pairs_eps[i][j] = sqrt(d_AtomList[i].epsilon * d_AtomList[j].epsilon);
            // s_ij = (s_i + s_j) / 2
            d_pairs_sig[i][j] = (d_AtomList[i].sigma + d_AtomList[j].sigma) / 2.0;
        }
    }
}
