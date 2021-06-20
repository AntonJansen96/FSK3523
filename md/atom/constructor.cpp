#include "atom.h"

Atom::Atom(size_t idx, std::string name, real mass, real epsilon, real sigma, std::vector<real> const &x)
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
