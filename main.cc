#include "md/md.h"

int main()
{
    // https://www.thermodynamik.tu-berlin.de/fileadmin/fg103/Publikationen/2017_nobleGases_LJ.pdf

    // std::vector<Atom> AtomList = {
    //     Atom(1, "ARG", 39.948, 0.998, 0.34, {1.6, 1.6, 1.5}),
    //     Atom(2, "NEO", 20.180, 0.282, 0.28, {2.2, 2.1, 1.5}),
    //     Atom(3, "ARG", 39.948, 0.998, 0.34, {2.4, 1.6, 1.5}),
    // };

    size_t idx = 1;
    std::vector<Atom> AtomList;
    for (double x = 1; x < 60; x += 5)
        for (double y = 1; y < 60; y += 5)
            for (double z = 1; z < 60; z += 5)
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
        500000,     // nsteps
        500,        // nstout
        0.0002,     // dt (ps)
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
