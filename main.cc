#include "atom/atom.h"
#include "md/md.h"

int main()
{
    std::vector<Atom> AtomList = {
        Atom(1, "ARG", 39.948, 0.998, 0.34, {0.1, 0.1, 0}),
        Atom(2, "ARG", 39.948, 0.998, 0.34, {0.5, 0.1, 0}),
        Atom(3, "ARG", 39.948, 0.998, 0.34, {0.9, 0.1, 0}),
    };

    MD simulator
    (
        10'000,     // nsteps
        100,        // nstout
        0.002,      // dt (ps)
        300,        // T (K)
        2,          // tauT
        1.2,        // LJcutoff (nm)
        false,      // useThermostat
        true,       // usePBC
        AtomList,   // AtomList
        {3, 3, 3}   // boxsize (nm)
    );
    
    simulator.run();
}
