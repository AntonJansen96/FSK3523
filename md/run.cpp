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
    Stopwatch timer1{"Generating velocities"};
    Stopwatch timer2{"MD-loop"};
    Stopwatch timer3{"Computing forces"};
    Stopwatch timer4{"Integration"};
    Stopwatch timer5{"Output"}; timer5.reset();

    // Pre-compute all the possible epsilon and sigma pairs.
    this->precomputepairs();

    // Generate velocities.
    timer1.start(); this->generate_velocities(); timer1.stop();

    // Output at step 0.
    this->writeFrame(0); 
    this->writeEnergies(0);

    // Do MD loop.
    timer2.start();

    for (size_t step = 1; step != d_nsteps + 1; ++step)
    {
        timer3.start();

        // Compute force.
        this->computeforces(step);
        
        timer3.stop(); timer4.start();

        // Update configuration.
        this->integrate();
        
        timer4.stop(); timer5.start();

        // Sample.
        if (step % d_nstout == 0)
        {
            this->writeFrame(step);
            this->writeEnergies(step);

            userUpdate(step, d_nsteps);
        }

        timer5.stop(); // Sample.
    }
    
    timer2.stop(); // MD-loop.

    printf("\n\n");
    timer1.time();    // Print generate velocities timing.
    timer2.time();    // Print MD-loop timing.
    timer3.time();    // Print computing forces timing.
    timer4.time();    // Print integration timing.
    timer5.time();    // Print output timing.
}
