/**
This header-only library contains simple templates for abstract containers to 
make life easier.
**/

#ifndef EASY_H
#define EASY_H

#include <numeric>      // for std::accumulate
#include <iostream>     // for std::cout
#include <math.h>       // for sqrt()

namespace easy {

template <typename Type>

// Returns sum of elements in abstract container.
typename Type::value_type sum(Type const &input)
{
    return std::accumulate
    (
        input.begin(), input.end(), typename Type::value_type()
    );
}

template <typename Type>

// Returns product of elements in abstract container. 
typename Type::value_type prod(Type const &input)
{
    typename Type::const_iterator iter, end(input.end());
    typename Type::value_type product = 1;

    for (iter = input.begin(); iter != end; ++iter)
        product *= *iter;

    return product;
}

template <typename Type>

// Returns mean of elements in abstract container.
double mean(Type const &input)
{
    return sum(input) / static_cast<double>(input.size());
}

template <typename Type>

// Returns the variance (sample) of abstract container.
double var(Type const &input)
{
    typename Type::const_iterator iter, end(input.end());
    double const mu = mean(input);

    double var = 0.0;
    for (iter = input.begin(); iter != end; ++iter)
        var += (*iter - mu) * (*iter - mu);

    return var / (input.size() - 1);
}

template <typename Type>

// Returns the standard deviation (sample) of abstract container.
double sdev(Type const &input)
{
    return std::sqrt(var(input));
}

template <typename Type>

// Return the largest element in the abstract container.
typename Type::value_type largest(Type const &input)
{
    typename Type::const_iterator iter, end(input.end());
    typename Type::value_type largest = *input.begin();

    for (iter = input.begin(); iter != end; ++iter)
        if (*iter > largest)
            largest = *iter;
    
    return largest;
}

template <typename Type>

// Return the smallest element in the abstract container.
typename Type::value_type smallest(Type const &input)
{
    typename Type::const_iterator iter, end(input.end());
    typename Type::value_type smallest = *input.begin();

    for (iter = input.begin(); iter != end; ++iter)
        if (*iter < smallest)
            smallest = *iter;
    
    return smallest;
}

template <typename Type>

// Print elements of abstract container to terminal.
void print(Type const &input)
{
    typename Type::const_iterator iter, end(input.end());

    for (iter = input.begin(); iter != end; ++iter)
        std::cout << *iter << ' ';

    std::cout << '\n';
}

} // Namespace.

#endif
