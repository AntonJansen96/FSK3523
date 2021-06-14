#include "md.h"

void MD::run()
{
    // Generate velocities.
    this->generate_velocities();

    // Output at step 0.
    this->writeFrame(0);
    this->writeEnergies(0);

    // Do MD loop.
    for (size_t step = 1; step != d_nsteps + 1; ++step)
    {
        // Compute force.
        this->get_accelerations(step);

        // Update configuration.
        this->integrate();

        // Sample.
        if (step % d_nstout == 0)
        {
            this->writeFrame(step);
            this->writeEnergies(step);
        }
    }
}
