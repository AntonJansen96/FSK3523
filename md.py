#!/bin/python3

# MODULES
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants.constants import Boltzmann

# CONSTANTS
m_argon = 39.948
k_b = 1.380649e-23

# FORCES
def lennard_jones_force(r, epsilon, sigma):
    repulsive  = 48 * epsilon * sigma**12 / r**13
    attractive = 24 * epsilon * sigma**6  / r**7

    return repulsive - attractive

def init_velocity(T, numParticles, mass):
    Boltzmann_factor = ((k_b * T) / (mass * 1.602e-19))**0.5
    
    velocities = numParticles * [0]
    for idx in range(0, numParticles):
        velocities[idx] = Boltzmann_factor * (np.random.rand() - 0.5)

    return velocities

def get_accelerations(positions):
    def reduce(forcegrid):
        reduced = len(forcegrid) * [0]

        for i in range(0, len(forcegrid)):
            for j in range(0, len(forcegrid)):
                reduced[i] += forcegrid[j][i]

        return reduced

    positions = list(positions)
    
    accel_x = [[0] * len(positions) for i in range(len(positions))]
        
    for i in range(0, len(positions) - 1):
        for j in range(i + 1, len(positions)):
            
            r_x = positions[j] - positions[i]
            
            rmag = (r_x**2)**0.5
            
            force_scalar = lennard_jones_force(rmag, 0.0103, 3.4)
            
            force_x = force_scalar * r_x / rmag
            
            accel_x[i][j] =   force_x / m_argon
            accel_x[j][i] = - force_x / m_argon

    accels = reduce(accel_x)
    
    return accels

def update_pos(x, v, a, dt):
    x_new = len(x) * [0]
    
    for i in range(0, len(x)):
        x_new[i] = x[i] + v[i] * dt + 0.5 * a[i] * dt * dt

    return x_new

def update_velo(v, a, a1, dt):
    v_new = len(v) * [0]

    for i in range(0, len(v)):
        v_new[i] = v[i] + 0.5 * (a[i] + a1[i]) * dt

    return v_new

def run_md(dt, number_of_steps, initial_temp, x, mass):
    positionList = [x]

    # Generate velocities and first acceleration    
    v = init_velocity(initial_temp, len(x), mass)
    a = get_accelerations(x)
    
    for i in range(number_of_steps):
        x = update_pos(x, v, a, dt)
        a1 = get_accelerations(x)
        v = update_velo(v, a, a1, dt)
        a = a1

        positionList.append(x)

    return positionList

# MAIN #########################################################################

# pos            dt  nsteps  T       coordinates      argonmass                                    
sim_pos = run_md(0.1, 10000, 300, [1, 5, 10], m_argon)

# PLOTTING #####################################################################

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
