#include "stopwatch.h"

// Move assignment.
Stopwatch &Stopwatch::operator=(Stopwatch &&tmp)
{
    if (not tmp.d_stopped)
        throw(std::logic_error("cannot move running Stopwatch"));
    
    swap(tmp);
    return *this;
}
