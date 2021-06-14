#ifndef MD_H
#define MD_H

#include "../atom/atom.h"
#include "constants.cc"
#include <vector>
#include <random>

class MD
{
    size_t d_nsteps;
    size_t d_nstout;
    double d_dt;
    double d_T;
    double d_tauT;
    double d_LJcutoff;
    bool d_useThermostat;
    bool d_usePBC;
    std::vector<Atom> d_AtomList;
    std::vector<double> d_boxsize;
    std::default_random_engine d_engine{1};
    double d_log_LJ_trun_energy = 0;
    double d_log_LJ_tail_energy = 0;

    public:
        // Constructor.
        MD
        (
            size_t nsteps, size_t nstout, double dt, double T, double tauT, 
            double LJcutoff, bool useThermostat, bool usePBC, 
            std::vector<Atom> const &AtomList, 
            std::vector<double> const &boxsize
        );

        // Run the molecular dynamics simulation.
        void run();
    
    private:
        // Generate initial velocities from Maxwell-Boltzmann distribution.
        void generate_velocities();
        
        // Compute the forces/accelerations acting on the particles.
        void get_accelerations(size_t step);
        
        // Update positions, velocities, and accelerations.
        void integrate();
        
        // Write trajectory frame to traj.pdb.
        void writeFrame(size_t step) const;
        
        // Write the energies to energy.log.
        void writeEnergies(size_t step) const;
};

#endif
