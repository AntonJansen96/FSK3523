#include "main.h"

// https://www.thermodynamik.tu-berlin.de/fileadmin/fg103/Publikationen/2017_nobleGases_LJ.pdf

// Generate a 3D lattice of particles for a given spacing and box size.
std::vector<Atom> genlattice(double spacing, std::vector<size_t> const &boxsize)
{
    size_t idx = 1;
    std::vector<Atom> Atoms;

    for (double x = 0.1; x < boxsize[0]; x += spacing)
    {
        for (double y = 0.1; y < boxsize[1]; y += spacing)
        {
            for (double z = 0.1; z < boxsize[2]; z += spacing)
            {
                if (idx % 2 == 0)
                    Atoms.push_back(Atom{idx, "ARG", 39.948, 0.998, 0.34, {x, y, z}});

                else
                    Atoms.push_back(Atom{idx, "NEO", 20.180, 0.282, 0.28, {x, y, z}});
                
                ++idx;
            }
        }
    }

    return Atoms;
}
