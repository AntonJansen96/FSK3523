#include "md.h"

// Idea: if we put the update of accelerations in the same loop as the
// reduced forces, we don't even need a vector to store the reduced
// values in, we can just use doubles. Note: this loop is thread safe.

void MD::reduceforces()
{
    #pragma omp parallel for
    for (size_t i = 0; i < d_AtomList.size(); ++i)
    {
        real Fx = 0;
        real Fy = 0;
        real Fz = 0;

        for (size_t j = 0; j < d_AtomList.size(); ++j)
        {
            Fx += d_Fgrid_x[j][i];
            Fy += d_Fgrid_y[j][i];
            Fz += d_Fgrid_z[j][i];
        }

        // Update the accelerations.
        d_AtomList[i].a_new[0] = Fx / d_AtomList[i].mass;
        d_AtomList[i].a_new[1] = Fy / d_AtomList[i].mass;
        d_AtomList[i].a_new[2] = Fz / d_AtomList[i].mass;
    }
}
