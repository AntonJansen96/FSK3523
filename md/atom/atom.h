#ifndef ATOM_H
#define ATOM_H

#include <vector>
#include <string>
#include "../../precision.h"

class Atom
{
    public:
        size_t      const idx;
        std::string const name;
        real        const mass;
        real        const epsilon;
        real        const sigma;
        std::vector<real> x;
        std::vector<real> v;
        std::vector<real> a;
        std::vector<real> a_new;
        
        // Constructor.
        Atom(size_t idx, std::string name, real mass, real epsilon, real sigma, std::vector<real> const &x);

        // Return magnitude of the velocity of the atom squared.
        real vmag_sq() const;
        
        // Return magnitude of the acceleration of the atom squared.
        real amag_sq() const;
        
        // Return magnitude of the acceleration_new of the atom squared.
        real a_new_mag_sq() const;
};

inline real Atom::vmag_sq() const
{
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

inline real Atom::amag_sq() const
{
    return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
}

inline real Atom::a_new_mag_sq() const
{
    return a_new[0] * a_new[0] + a_new[1] * a_new[1] + a_new[2] * a_new[2];
}

#endif
