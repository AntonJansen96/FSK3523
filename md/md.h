#ifndef MD_H
#define MD_H

// #define DEBUG
#ifdef DEBUG
#include "easy.h"
#endif

#include "atom/atom.h"
#include <vector>
#include <random>

typedef std::vector<std::vector<double>> grid;

class MD
{
    // Molecular dynamics parameters.
    size_t const d_nsteps;
    size_t const d_nstout;
    double const d_dt;
    double const d_T;
    double const d_tauT;
    double const d_LJcutoff;
    bool   const d_useLJ;
    bool   const d_useThermostat;
    bool   const d_usePBC;

    // The structure and topology.
    std::vector<Atom> d_AtomList;
    std::vector<double> const d_boxsize;
    std::default_random_engine d_engine{1};

    // Stores pre-computed combinations of epsilon and sigma.
    grid d_pairs_eps;
    grid d_pairs_sig;

    // For energy logging.
    double d_log_LJ_energy = 0;
    double d_log_thermo_energy = 0;

    // For speedup of computeforces().
    grid const d_emptyFgrid;
    grid d_Fgrid_x, d_Fgrid_y, d_Fgrid_z;

    public:
        // Constructor.
        MD
        (
            size_t nsteps, size_t nstout, double dt, double T, double tauT, 
            double LJcutoff, bool useLJ, bool useThermostat, bool usePBC, 
            std::vector<Atom> const &AtomList, 
            std::vector<double> const &boxsize
        );

        // This object is not meant to be copied or moved.
        MD(MD const &other) = delete;
        MD(MD &&temp) = delete;

        // Run the molecular dynamics simulation.
        void run();

    private:
        // Generate initial velocities from Maxwell-Boltzmann distribution.
        void generate_velocities();
        
        // Compute the forces/accelerations acting on the particles.
        void computeforces(size_t step);

        // Pre-compute combinations of epsilon and sigma.
        void precomputepairs();

        // Compute Lennard-Jones tail-correction to energy.
        double tailcorrection() const;

        // Reduces the forces from the force grids and updates the accelerations.
        void reduceforces();

        // Velocity-Verlet integrator. Update positions, velocities, and accelerations.
        void integrate();
        
        // Write trajectory frame to traj.pdb.
        void writeFrame(size_t step) const;
        
        // Write the energies to energy.log.
        void writeEnergies(size_t step) const;
};

#endif
