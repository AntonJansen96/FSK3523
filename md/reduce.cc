#include "md.h"

std::vector<double> MD::reduce(std::vector<std::vector<double>> const &forcegrid)
{
    std::vector<double> reduced(0, forcegrid.size());

    for (size_t i = 0; i < forcegrid.size(); ++i)
    {
        for (size_t j = 0; j < forcegrid.size(); ++j)
        {
            reduced[i] += forcegrid[j][i];
        }
    }

    return reduced;
}
