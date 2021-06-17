#ifndef ATOM_H
#define ATOM_H

#include <vector>
#include <string>

class Atom
{
    public:
        size_t      const idx;
        std::string const name;
        double      const mass;
        double      const epsilon;
        double      const sigma;
        std::vector<double> x;
        std::vector<double> v;
        std::vector<double> a;
        std::vector<double> a_new;
        
        // Constructor.
        Atom(size_t idx, std::string name, double mass, double epsilon, double sigma, std::vector<double> const &x);

        // Return magnitude of the velocity of the atom squared.
        double vmag_sq() const;
        
        // Return magnitude of the acceleration of the atom squared.
        double amag_sq() const;
        
        // Return magnitude of the acceleration_new of the atom squared.
        double a_new_mag_sq() const;
};

inline double Atom::vmag_sq() const
{
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

inline double Atom::amag_sq() const
{
    return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
}

inline double Atom::a_new_mag_sq() const
{
    return a_new[0] * a_new[0] + a_new[1] * a_new[1] + a_new[2] * a_new[2];
}

#endif
