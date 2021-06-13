#include "easy.h"
#include <iostream>
#include <vector>
#include <string>

namespace constants
{
    double const N   = 6.02214129e23;    // mol^-1
    double const kB  = 1.380649e-23;     // J/K
    double const amu = 1.660539066e-27;  // kg
    double const pi  = 3.141592654;
}

class Atom
{
    public:
        // Data members
        size_t      idx;
        std::string name;
        double      mass;
        double      epsilon;
        double      sigma;
        
        std::vector<double> x;
        std::vector<double> v;
        std::vector<double> a;
        std::vector<double> a_new;

        // Constructor
        Atom
        (   
            size_t idx, std::string name, double mass, double epsilon, 
            double sigma, std::vector<double> const &x
        )
        :
            idx(idx), mass(mass), epsilon(epsilon), sigma(sigma), x(x),
            v({0, 0, 0}), a({0, 0, 0}), a_new({0, 0, 0})
        {}
};

class MD
{
    // Data members (private).
    size_t nsteps;
    size_t nstout;
    double dt;
    double T;
    double tauT;
    double LJcutoff;
    bool useThermostat;
    bool usePBC;
    std::vector<Atom> AtomList;
    std::vector<double> boxsize;

    public:
        void generate_velocities();
        void get_accelerations();
        void integrate();
        void writeFrame();
        void writeEnergies();
    
    private:
        std::vector<double> reduce();
};

int main()
{
    std::vector<Atom> AtomList = {
        Atom(1, "ARG", 39.948, 0.998, 0.34, {1, 1, 0}),
        Atom(2, "ARG", 39.948, 0.998, 0.34, {5, 1, 0}),
        Atom(3, "ARG", 39.948, 0.998, 0.34, {9, 1, 0}),
    };

    std::cout << AtomList[0].mass << '\n';
    std::cout << constants::amu << '\n';
}
