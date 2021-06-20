#include "md.h"

// Constructor.
MD::MD
(
    size_t nsteps, size_t nstout, double dt, double T, double tauT, 
    double LJcutoff, bool useLJ, bool useThermostat, bool usePBC, 
    std::vector<Atom> const &AtomList, 
    std::vector<double> const &boxsize
)
:
    d_nsteps(nsteps),
    d_nstout(nstout),
    d_dt(dt),
    d_T(T),
    d_tauT(tauT),
    d_LJcutoff(LJcutoff),
    d_useLJ(useLJ),
    d_useThermostat(useThermostat),
    d_usePBC(usePBC),
    d_AtomList(AtomList),
    d_boxsize(boxsize),
    d_emptyFgrid(grid(d_AtomList.size(), std::vector<double>(d_AtomList.size(), 0)))
{
    printf("\nSuccesfully initialized MD object (%zu atoms):\n", d_AtomList.size());
    
    printf("nsteps        = %zu\n"  , d_nsteps);
    printf("nstout        = %zu\n"  , d_nstout);
    printf("dt            = %.4f ps\n", d_dt);
    printf("T             = %.1f K\n" , d_T);
    printf("tauT          = %.2f\n"   , d_tauT);
    printf("LJcutoff      = %.2f nm\n", d_LJcutoff);
    printf("useLJ         = %d\n"   , d_useLJ);
    printf("useThermostat = %d\n"   , d_useThermostat);
    printf("usePBC        = %d\n"   , d_usePBC);
    printf("box           = [%.2f, %.2f, %.2f]\n\n", d_boxsize[0], d_boxsize[1], d_boxsize[2]);
}
