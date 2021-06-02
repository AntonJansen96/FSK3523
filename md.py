#!/bin/python3

import numpy as np
import matplotlib.pyplot as plt

class constants:
    k_b = 1.380649e-23    # J/K
    e   = 1.602176634e-19 # C

def lennard_jones_force(r, epsilon, sigma):
    repulsive  = 48 * epsilon * sigma**12 / r**13
    attractive = 24 * epsilon * sigma**6  / r**7

    return repulsive - attractive

def init_velocity(T, numParticles, mass):
    Boltzmann_factor = ((constants.k_b * T) / (mass * constants.e))**0.5
    
    velocities = numParticles * [0]
    for idx in range(0, numParticles):
        velocities[idx] = Boltzmann_factor * (np.random.rand() - 0.5)

    return velocities

def get_accelerations(positions, mass):
    def reduce(forcegrid):
        reduced = len(forcegrid) * [0]

        for i in range(0, len(forcegrid)):
            for j in range(0, len(forcegrid)):
                reduced[i] += forcegrid[j][i]

        return reduced
    
    accel_x = [[0] * len(positions) for i in range(len(positions))]
        
    for i in range(0, len(positions) - 1):
        for j in range(i + 1, len(positions)):
            
            r_x = positions[j] - positions[i]
            
            rmag = (r_x**2)**0.5
            
            force_scalar = lennard_jones_force(rmag, 0.0103, 3.4)
            
            force_x = force_scalar * r_x / rmag
            
            accel_x[i][j] =   force_x / mass
            accel_x[j][i] = - force_x / mass

    accels = reduce(accel_x)
    
    return accels

# Velocity-Verlet integrator.
def integrator(x, v, a, a_new, dt):
    x_new = len(x) * [0]
    v_new = len(x) * [0]
    
    for i in range(0, len(x)):
        x_new[i] = x[i] + v[i] * dt + 0.5 * a[i] * dt * dt
        v_new[i] = v[i] + 0.5 * (a[i] + a_new[i]) * dt

    return x_new, v_new, a_new

# Run MD simulation.
def run_md(dt, nsteps, T, x, mass):
    positionList = [x]

    # 1 INPUT INITIAL CONDITIONS
    
    # Initial velocity from Maxwell-Boltzmann distribution.
    v = init_velocity(T=T, numParticles=len(x), mass=mass)
    
    # Initial acceleration is zero.
    a = len(x) * [0]
    
    for _ in range(nsteps):

        # 2 COMPUTE FORCE
        a_new = get_accelerations(x, mass)

        # 3 UPDATE CONFIGURATION
        x, v, a = integrator(x, v, a, a_new, dt)

        # 4 OUTPUT STEP
        positionList.append(x)

    return positionList

# MAIN #########################################################################

sim_pos = run_md(dt=0.1, nsteps=10000, T=300, x=[1, 5, 10], mass=39.948)

atom1 = len(sim_pos) * [0]
atom2 = len(sim_pos) * [0]
atom3 = len(sim_pos) * [0]

for i in range(0, len(sim_pos)):
    atom1[i] = sim_pos[i][0]
    atom2[i] = sim_pos[i][1]
    atom3[i] = sim_pos[i][2]

plt.plot(atom1, '.', label='atom 1')
plt.plot(atom2, '.', label='atom 2')
plt.plot(atom3, '.', label='atom 3')

plt.xlabel(r'Step')
plt.ylabel(r'$x$-Position/Å')
plt.legend(frameon=False)
plt.show()
