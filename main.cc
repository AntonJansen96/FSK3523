#include "md/md.h"

int main()
{
    // https://www.thermodynamik.tu-berlin.de/fileadmin/fg103/Publikationen/2017_nobleGases_LJ.pdf

    // std::vector<Atom> AtomList = {
    //     Atom(1, "ARG", 39.948, 0.998, 0.34, {0.1, 0.1, 0}),
    //     Atom(2, "ARG", 39.948, 0.998, 0.34, {0.5, 0.1, 0}),
    //     Atom(3, "ARG", 39.948, 0.998, 0.34, {0.9, 0.1, 0}),
    // };

    size_t idx = 1;
    std::vector<Atom> AtomList;
    for (double x = 1; x < 30; x += 8)
        for (double y = 1; y < 30; y += 8)
            for (double z = 1; z < 30; z += 8)
            {
                if (idx % 2 == 0)
                {
                AtomList.push_back
                (
                    Atom(idx, "ARG", 39.948, 0.998, 0.34, {x/10, y/10, z/10})
                );
                }
                else
                {
                    AtomList.push_back
                    (
                    Atom(idx, "NEO", 20.180, 0.282, 0.28, {x/10, y/10, z/10})
                    );
                }

                ++idx;
            }

    printf("we have %zu atoms...\n\n", AtomList.size());

    MD simulator
    (
        10'000,     // nsteps
        100,        // nstout
        0.002,      // dt (ps)
        300,        // T (K)
        1,          // tauT
        1.2,        // LJcutoff (nm)
        true,       // useLJ
        true,       // useThermostat
        true,       // usePBC
        AtomList,   // AtomList
        {6, 6, 6}   // boxsize (nm)
    );
    
    simulator.run();
}
