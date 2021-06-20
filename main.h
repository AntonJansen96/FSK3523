#include "md/md.h"
#include <omp.h>

// Generate a 3D lattice of particles for a given spacing and box size.
std::vector<Atom> genlattice(double spacing, std::vector<size_t> const &boxsize);

// Lattice for testing.
std::vector<Atom> const testLattice = 
{
    Atom(1, "ARG", 39.948, 0.998, 0.34, {1.6, 1.6, 1.5}),
    Atom(2, "NEO", 20.180, 0.282, 0.28, {2.2, 2.1, 1.5}),
    Atom(3, "ARG", 39.948, 0.998, 0.34, {2.4, 1.6, 1.5}),
};
