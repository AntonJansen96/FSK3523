#include "md.h"
#include "stopwatch/stopwatch.h"

namespace {
void userUpdate(size_t step, size_t d_nsteps)
{
    printf
    (
        "\rstep %zu, %.1f%% (%zu)", 
        step, 
        100 * step/static_cast<real>(d_nsteps), d_nsteps
    );
    
    fflush(stdout);
}
} // namespace

void MD::run()
{
    // Initialize Stopwatch objects for timing.
    Stopwatch timer1{"MD-loop"};    
    Stopwatch timer2{"Computing forces"};
    Stopwatch timer3{"Integration"};
    Stopwatch timer4{"Output"}; timer4.reset();

    // Pre-compute all the possible epsilon and sigma pairs.
    this->precomputepairs();

    // Generate velocities.
    this->generate_velocities();

    // Do MD-loop.
    timer1.start();
    for (size_t step = 0; step != d_nsteps + 1; ++step)
    {
        // Compute forces.
        timer2.start(); this->computeforces(step); timer2.stop();

        // Sample.
        timer4.start();
        if (step % d_nstout == 0)
        {
            this->writeFrame(step);
            this->writeEnergies(step);

            userUpdate(step, d_nsteps);
        } timer4.stop();

        // Update velocities and positions.
        timer3.start(); this->integrate(); timer3.stop();
    } timer1.stop();

    printf("\n\n");
    timer1.time();    // Print MD-loop timing.
    timer2.time();    // Print computing forces timing.
    timer3.time();    // Print integration timing.
    timer4.time();    // Print output timing.
}
