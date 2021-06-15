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

        Atom(size_t idx, std::string name, double mass, double epsilon, double sigma, std::vector<double> const &x);
};

#endif
