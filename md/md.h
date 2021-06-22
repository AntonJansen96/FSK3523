#ifndef MD_H
#define MD_H

#include "../precision.h"
#include "atom/atom.h"
#include <random>

typedef std::vector<std::vector<real>> grid;

class MD
{
    // Molecular dynamics parameters.
    size_t const d_nsteps;
    size_t const d_nstout;
    real   const d_dt;
    real   const d_T;
    real   const d_tauT;
    real   const d_LJcutoff;
    bool   const d_useLJ;
    bool   const d_useThermostat;
    bool   const d_usePBC;

    // The structure and topology.
    std::vector<Atom> d_AtomList;
    std::vector<real> const d_boxsize;
    std::default_random_engine d_engine{1};

    // Stores pre-computed combinations of epsilon and sigma.
    grid d_pairs_eps;
    grid d_pairs_sig;

    // For energy logging.
    real d_log_LJ_energy = 0;
    real d_log_thermo_energy = 0;

    // For speedup of computeforces().
    grid const d_emptyFgrid;
    grid d_Fgrid_x, d_Fgrid_y, d_Fgrid_z;
    real const d_LJcutoff_6;

    public:
        // Constructor.
        MD
        (
            size_t nsteps, size_t nstout, real dt, real T, real tauT, 
            real LJcutoff, bool useLJ, bool useThermostat, bool usePBC, 
            std::vector<Atom> const &AtomList, 
            std::vector<real> const &boxsize
        );

        // This object is not meant to be copied or moved.
        MD(MD const &other) = delete;
        MD(MD &&temp) = delete;

        // Run the molecular dynamics simulation.
        void run();

    private:
        // Pre-compute combinations of epsilon and sigma.
        void precomputepairs();

        // Generate initial velocities from Maxwell-Boltzmann distribution.
        void generate_velocities();

        // Velocity-Verlet step one. Update positions and apply PBC.
        void velocityVerlet1();

        // Compute the forces acting on the particles.
        void computeforces(size_t step);

        // Reduces the forces from the force grids and updates the accelerations.
        void reduceforces();

        // Velocity-Verlet step two. Update velocities and accelerations.
        void velocityVerlet2();

        // Write trajectory frame to traj.pdb.
        void writeFrame(size_t step) const;
        
        // Write the energies to energy.log.
        void writeEnergies(size_t step) const;

        // Compute Lennard-Jones tail-correction to energy.
        real tailcorrection() const;
};

#endif
