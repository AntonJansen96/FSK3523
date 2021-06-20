#include "main.h"

int main()
{
    omp_set_num_threads(8); // Multithreading.

    std::vector<Atom> const AtomList = genlattice(0.5, {3, 3, 3});

    MD
    (
        50'000,     // nsteps
        100,        // nstout
        0.0002,     // dt (ps)
        300,        // T (K)
        1,          // tauT
        1.2,        // LJcutoff (nm)
        true,       // useLJ
        true,       // useThermostat
        true,       // usePBC
        AtomList,   // AtomList
        {3, 3, 3}   // boxsize (nm)
    ).run();
}
