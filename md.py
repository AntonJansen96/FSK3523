#!/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import time
import load
import math

# BASE UNITS
# length    nm      1e-9 m
# mass      amu     1.66e-27 kg
# time      ps      1e-12 s
# charge    e       1.602e-19 C
# temp      K

# DERIVED UNITS
# velocity  length/time     nm/ps       1e3 m/s
# energy                    amu nm^2 / ps^2 = kJ/mol

class constants:
    N   = 6.02214129e23   # mol^-1
    kB  = 1.380649e-23    # J/K
    R   = N * kB          # J/mol/K  
    e0  = 8.854187813e-12 # Farad/m
    e   = 1.602176634e-19 # C
    amu = 1.660539066e-27 # kg

class Atom:
    def __init__(self, idx, name, x, mass, charge, lj_epsilon, lj_sigma):
        self.idx        = idx
        self.name       = name
        self.mass       = mass
        self.charge     = charge
        self.lj_epsilon = lj_epsilon
        self.lj_sigma   = lj_sigma
        self.x          = x
        self.v          = [0, 0, 0]
        self.a          = [0, 0, 0]
        self.a_new      = [0, 0, 0]

def generate_velocities(AtomList, T):
    # We go from m/s to nm/ps so we have a factor 0.001. Furthermore, we have 
    # a sqrt(2) from 0.5mv^2 and 1/sqrt(3) from vector calculus.
    factor = 0.001 * math.sqrt(2.0 / 3.0)
    
    for Atom in AtomList:
        sigma = math.sqrt((constants.kB * T) / (Atom.mass * constants.amu))

        Atom.v[0] = factor * np.random.normal(scale=sigma)
        Atom.v[1] = factor * np.random.normal(scale=sigma)
        Atom.v[2] = factor * np.random.normal(scale=sigma)

def get_accelerations(AtomList, boxsize):
    def reduce(forcegrid):
        reduced = len(forcegrid) * [0]

        for i in range(0, len(forcegrid)):
            for j in range(0, len(forcegrid)):
                reduced[i] += forcegrid[j][i]

        return reduced

    accel_x = [[0] * len(AtomList) for i in range(len(AtomList))]
    accel_y = [[0] * len(AtomList) for i in range(len(AtomList))]
    accel_z = [[0] * len(AtomList) for i in range(len(AtomList))]

    # Loop over all combinations.
    for i in range(0, len(AtomList) - 1):
        for j in range(i + 1, len(AtomList)):

            # Get the distance in each dimension.
            r_x = AtomList[j].x[0] - AtomList[i].x[0]
            r_y = AtomList[j].x[1] - AtomList[i].x[1]
            r_z = AtomList[j].x[2] - AtomList[i].x[2]

            # Apply minimum image convention.
            if (r_x >     0.5 * boxsize[0]): r_x -= boxsize[0]
            elif (r_x <= -0.5 * boxsize[0]): r_x += boxsize[0]

            if (r_y >     0.5 * boxsize[1]): r_y -= boxsize[1]
            elif (r_y <= -0.5 * boxsize[1]): r_y += boxsize[1]

            if (r_z >     0.5 * boxsize[2]): r_z -= boxsize[2]
            elif (r_z <= -0.5 * boxsize[2]): r_z += boxsize[2]

            # Compute distance.
            rmag_sq = r_x * r_x + r_y * r_y + r_z * r_z

            # Compute Lennard-Jones force.
            if (rmag_sq < 1.44): # Cut-off = 1.2nm, hardcoded.
                rmag = math.sqrt(rmag_sq)
                
                # Combination rule (arithmetic mean).
                epsilon = 0.5 * (AtomList[i].lj_epsilon + AtomList[j].lj_epsilon) # in kJ/mol (= natural MD unit)
                sigma   = 0.5 * (AtomList[i].lj_sigma + AtomList[j].lj_sigma)     # in nm

                # Compute LJ terms.
                factor       = (sigma * sigma * sigma * sigma * sigma * sigma) / (rmag_sq * rmag_sq * rmag_sq)
                attractive   = factor / rmag
                repulsive    = attractive * factor
                force_scalar = 24 * epsilon * (2 * repulsive - attractive)

                # Decompose force.
                force_x = force_scalar * r_x / rmag
                force_y = force_scalar * r_y / rmag
                force_z = force_scalar * r_z / rmag

                # Fill the force grid (are the indices of the masses correct?).
                accel_x[i][j] =   force_x / AtomList[i].mass
                accel_x[j][i] = - force_x / AtomList[j].mass

                accel_y[i][j] =   force_y / AtomList[i].mass
                accel_y[j][i] = - force_y / AtomList[j].mass

                accel_z[i][j] =   force_z / AtomList[i].mass
                accel_z[j][i] = - force_z / AtomList[j].mass
            
            else:
                # Fill the force grid.
                accel_x[i][j] = 0; accel_x[j][i] = 0
                accel_y[i][j] = 0; accel_y[j][i] = 0
                accel_z[i][j] = 0; accel_z[j][i] = 0

    # Reduce forces.
    reduced_x = reduce(accel_x)
    reduced_y = reduce(accel_y)
    reduced_z = reduce(accel_z)

    # Update.
    for i in range(0, len(AtomList)):
        AtomList[i].a_new = [reduced_x[i], reduced_y[i], reduced_z[i]]

def instaTemp(AtomList):
    # instantaneous temperature = sum(mv^2) / 2Nk_B
    # v is in nm/ps = 1e3 m/s --> v^2 = (nm/ps)^2 = 1e6 (m/s)^2
    # m is in amu = 1.66e-27 kg
    # factor = (10**6 * constants.amu) / (2 * constants.kB)
    # factor = 60.13618
    temp = 0
    for Atom in AtomList:
        vmag_sq = Atom.v[0] * Atom.v[0] + Atom.v[1] * Atom.v[1] + Atom.v[2] * Atom.v[2]
        temp += Atom.mass * vmag_sq

    return 60.13618 * temp / len(AtomList)

def kinetic_energy(AtomList):
    Ekin = 0
    for Atom in AtomList:
        vmag_sq = Atom.v[0] * Atom.v[0] + Atom.v[1] * Atom.v[1] + Atom.v[2] * Atom.v[2]
        # Use E = 1/2mv^2 and convert to J.
        Ekin += 0.5 * (constants.amu * Atom.mass) * (10**6 * vmag_sq)
    
    # convert J to kJ/mol.
    return 0.001 * Ekin * constants.N

def integrate(AtomList, dt, boxsize, T, tau_t):
    for Atom in AtomList:
        # Update positions.
        Atom.x[0] += (Atom.v[0] + 0.5 * Atom.a[0] * dt) * dt
        Atom.x[1] += (Atom.v[1] + 0.5 * Atom.a[1] * dt) * dt
        Atom.x[2] += (Atom.v[2] + 0.5 * Atom.a[2] * dt) * dt

        # Apply periodic boundary condition.
        if (Atom.x[0] > boxsize[0]): Atom.x[0] -= boxsize[0]
        elif (Atom.x[0] < 0       ): Atom.x[0] += boxsize[0]

        if (Atom.x[1] > boxsize[1]): Atom.x[1] -= boxsize[1]
        elif (Atom.x[1] < 0       ): Atom.x[1] += boxsize[1]

        if (Atom.x[2] > boxsize[2]): Atom.x[2] -= boxsize[2]
        elif (Atom.x[2] < 0       ): Atom.x[2] += boxsize[2]

        # Update velocities.
        if (np.random.rand() < tau_t * dt):
            # By applying Andersen thermostat:
            factor = 0.001 * math.sqrt(2.0 / 3.0)
            
            sigma = math.sqrt((constants.kB * T) / (Atom.mass * constants.amu))

            Atom.v[0] = factor * np.random.normal(scale=sigma)
            Atom.v[1] = factor * np.random.normal(scale=sigma)
            Atom.v[2] = factor * np.random.normal(scale=sigma)
        else:
            # By updating the velocities like usual:
            Atom.v[0] += 0.5 * (Atom.a[0] + Atom.a_new[0]) * dt
            Atom.v[1] += 0.5 * (Atom.a[1] + Atom.a_new[1]) * dt
            Atom.v[2] += 0.5 * (Atom.a[2] + Atom.a_new[2]) * dt

        # Update accelerations.
        Atom.a[0] = Atom.a_new[0]
        Atom.a[1] = Atom.a_new[1]
        Atom.a[2] = Atom.a_new[2]

# Write coordinates to a .pdb trajectory file
def writeFrame(step, AtomList, boxsize):
    with open("traj.pdb", "a+") as file:

        file.write("MODEL        {}\n".format(step))
        file.write("CRYST1  {:7.3f}  {:7.3f}  {:7.3f}  90.00  90.00  90.00 P 1           1\n".format(10 * boxsize[0], 10 * boxsize[1], 10 * boxsize[2]))

        for Atom in AtomList:
            file.write("{:6s}{:5d} {:^4s}{:1s}{:4s}{:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}\n".format('ATOM', Atom.idx, Atom.name, '', Atom.name, '', Atom.idx, '', 10 * Atom.x[0], 10 * Atom.x[1], 10 * Atom.x[2]))

        file.write('TER\nENDMDL\n')

# Function for logging various (ensemble) parameters in the system/
def writeEnergies(step, AtomList):
    with open("energy.log", "a+") as file:
        file.write("{} {} {}\n".format(step, instaTemp(AtomList), kinetic_energy(AtomList)))

# Run MD simulation.
def run_md(AtomList, dt, nsteps, T, boxsize):
    # 1 INPUT INITIAL CONDITIONS

        # Initialize list of positions for old output function.
    positionList = [[Atom.x[0] for Atom in AtomList]]

        # Initial velocity from Maxwell-Boltzmann distribution.
    generate_velocities(AtomList, T)

        # Update user.
    print("initial positions:  {:.4f}, {:.4f}, {:.4f}".format(AtomList[0].x[0], AtomList[1].x[0], AtomList[2].x[0]))
    print("initial velocities: {:.4f}, {:.4f}, {:.4f}".format(AtomList[0].v[0], AtomList[1].v[0], AtomList[2].v[0]))
    writeEnergies(0, AtomList) # Write the temperature at t = 0.

    time_forces    = 0
    time_integrate = 0
    time_output    = 0

    for step in range(nsteps):
        # 2 COMPUTE FORCE
        tic = time.perf_counter()
        get_accelerations(AtomList, boxsize)
        time_forces += time.perf_counter() - tic

        # 3 UPDATE CONFIGURATION
        tic = time.perf_counter()
        integrate(AtomList, dt, boxsize, T, tau_t=2)
        time_integrate += time.perf_counter() - tic

        # 4 OUTPUT STEP
        tic = time.perf_counter()
        if (step % 100 == 0):
            writeFrame(step + 1, AtomList, boxsize)
            writeEnergies(step + 1, AtomList)
            positionList.append([Atom.x[0] for Atom in AtomList])
        time_output += time.perf_counter() - tic

    print("time_forces = {:.4f}".format(time_forces))
    print("time_integrate = {:.4f}".format(time_integrate))
    print("time_output = {:.4f}".format(time_output))

    return positionList

# MAIN #########################################################################

np.random.seed(1) # for testing

# Argon has mass of 39.948 amu, an LJ-sigma of 0.34 nm, and an LJ-epsilon/kb = 
# 120 --> LJ-epsilon = 0.998 kJ/mol.

AtomList = [
    Atom(1, 'ARG', [1.0, 1, 0], 39.948, 0, 0.998, 0.34),
    Atom(2, 'ARG', [5.0, 1, 0], 39.948, 0, 0.998, 0.34),
    Atom(3, 'ARG', [9.0, 1, 0], 39.948, 0, 0.998, 0.34)
    ]

# AtomList = []
# idx = 1
# for x in range(1, 30, 8):
#     for y in range(1, 30, 8):
#         for z in range(1, 30, 8):
#             AtomList.append(Atom(idx, 'ARG', [x/10., y/10., z/10.], 39.948, 0, 0.998, 0.34))
#             idx += 1

print("we have {} atoms".format(len(AtomList)))

tic = time.perf_counter()
sim_pos = run_md(AtomList=AtomList, dt=0.002, nsteps=10000, T=300, boxsize=[3.0, 3.0, 3.0])
toc = time.perf_counter()
print("time_total = {:.4f}".format(toc - tic))

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

x = load.Col('energy.log', 1)
T = load.Col('energy.log', 2)

plt.plot(x, T, label="mean = {:.4f}".format(sum(T)/len(T)))
plt.legend()
plt.show()
