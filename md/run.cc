#include "md.h"
#include "stopwatch/stopwatch.h"

void MD::run()
{
    // Initialize Stopwatch objects for timing.
    Stopwatch timer1{"Generating velocities"};
    Stopwatch timer2{"MD-loop"};
    Stopwatch timer3{"Computing forces"};
    Stopwatch timer4{"Integration"};
    Stopwatch timer5{"Output"}; timer5.reset();

    // Generate velocities.
    timer1.start(); this->generate_velocities(); timer1.stop();

    // Output at step 0.
    this->writeFrame(0); 
    this->writeEnergies(0);

    // Do MD loop.
    timer2.start();
    for (size_t step = 1; step != d_nsteps + 1; ++step)
    {
        // Compute force.
        timer3.start(); this->get_accelerations(step); timer3.stop();

        // Update configuration.
        timer4.start(); this->integrate(); timer4.stop();

        // Sample.
        timer5.start();
        if (step % d_nstout == 0)
        {
            this->writeFrame(step);
            this->writeEnergies(step);
        } timer5.stop();
    } timer2.stop();

    timer1.time();    // Print generate velocities timing.
    timer2.time();    // Print MD-loop timing.
    timer3.time();    // Print computing forces timing.
    timer4.time();    // Print integration timing.
    timer5.time();    // Print output timing.
}
