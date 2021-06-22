#include "md.h"

void MD::velocityVerlet1()
{
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
                else if (atom.x[i] < 0.0)
                    atom.x[i] += d_boxsize[i];
            }
        }
    }
}
