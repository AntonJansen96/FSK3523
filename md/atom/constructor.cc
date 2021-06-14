#include "atom.h"

Atom::Atom(size_t idx, std::string name, double mass, double epsilon, double sigma, std::vector<double> const &x)
:
    idx(idx),
    name(name),
    mass(mass),
    epsilon(epsilon),
    sigma(sigma),
    x(x),
    v({0, 0, 0}),
    a({0, 0, 0}),
    a_new({0, 0, 0})
{}
