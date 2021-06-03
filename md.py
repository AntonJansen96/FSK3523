#!/bin/python3

import numpy as np
import matplotlib.pyplot as plt

class constants:
    N  = 6.02214129e23   # mol^-1
    kB = 1.380649e-23    # J/K
    R  = N * kB          # J/mol/K  
    e0 = 8.854187813e-12 # Farad/m
    e  = 1.602176634e-19 # C
    u  = 1.660539066e-27 # kg

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

def LJ_force(r, epsilon, sigma):
    LJ_cutoff = 12 # Angstrom = 1.2nm, hardcoded.
    
    if (r < LJ_cutoff):
        repulsive  = 48 * epsilon * sigma**12 / r**13
        attractive = 24 * epsilon * sigma**6  / r**7

        return repulsive - attractive

    return 0

def generate_velocities(AtomList, T):
    # Not sure if this is exactly correct.
    # Maybe some sqrt(3/2) is still missing?
    # See https://manual.gromacs.org/current/reference-manual/algorithms/molecular-dynamics.html#coordinates-and-velocities
    for Atom in AtomList:
        BoltzmannFactor = constants.kB * T
        BoltzmannFactor = BoltzmannFactor / (Atom.mass * constants.e)
        BoltzmannFactor = BoltzmannFactor**0.5

        Atom.v[0] = BoltzmannFactor * (np.random.rand() - 0.5)
        Atom.v[1] = BoltzmannFactor * (np.random.rand() - 0.5)
        Atom.v[2] = BoltzmannFactor * (np.random.rand() - 0.5)

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
        
    for i in range(0, len(AtomList) - 1):
        for j in range(i + 1, len(AtomList)):
            
            r_x = AtomList[j].x[0] - AtomList[i].x[0]
            r_y = AtomList[j].x[1] - AtomList[i].x[1]
            r_z = AtomList[j].x[2] - AtomList[i].x[2]

            # Periodic boundary condition
            r_x = r_x - np.round(r_x / boxsize[0]) * boxsize[0]
            r_y = r_y - np.round(r_y / boxsize[1]) * boxsize[1]
            r_z = r_z - np.round(r_z / boxsize[2]) * boxsize[2]

            # if r_x >   0.5 * boxsize[0]: r_x = r_x - boxsize[0]
            # if r_x <= -0.5 * boxsize[0]: r_x = r_x + boxsize[0]
            
            # if r_y >   0.5 * boxsize[1]: r_y = r_y - boxsize[1]
            # if r_y <= -0.5 * boxsize[1]: r_y = r_y + boxsize[1]
            
            # if r_z >   0.5 * boxsize[2]: r_z = r_z - boxsize[2]
            # if r_z <= -0.5 * boxsize[2]: r_z = r_z + boxsize[2]

            rmag = (r_x**2 + r_y**2 + r_z**2)**0.5
            
            force_scalar = LJ_force(rmag, AtomList[i].lj_epsilon, AtomList[i].lj_sigma)

            force_x = force_scalar * r_x / rmag
            force_y = force_scalar * r_y / rmag
            force_z = force_scalar * r_z / rmag
            
            accel_x[i][j] =   force_x / AtomList[i].mass
            accel_x[j][i] = - force_x / AtomList[j].mass

            accel_y[i][j] =   force_y / AtomList[i].mass
            accel_y[j][i] = - force_y / AtomList[j].mass

            accel_z[i][j] =   force_z / AtomList[i].mass
            accel_z[j][i] = - force_z / AtomList[j].mass

    reduced_x = reduce(accel_x)
    reduced_y = reduce(accel_y)
    reduced_z = reduce(accel_z)
    
    for i in range(0, len(AtomList)):
        AtomList[i].a_new[0] = reduced_x[i]
        AtomList[i].a_new[1] = reduced_y[i]
        AtomList[i].a_new[2] = reduced_z[i]

# Velocity-Verlet integrator.
def integrate(AtomList, dt, boxsize):
    for Atom in AtomList:
        Atom.x[0] += Atom.v[0] * dt + 0.5 * Atom.a[0] * dt * dt
        Atom.x[1] += Atom.v[1] * dt + 0.5 * Atom.a[1] * dt * dt
        Atom.x[2] += Atom.v[2] * dt + 0.5 * Atom.a[2] * dt * dt

        # Periodic boundary condition
        # if Atom.x[0] < -0.5 * boxsize[0]: Atom.x[0] = Atom.x[0] + boxsize[0]
        # if Atom.x[0] >= 0.5 * boxsize[0]: Atom.x[0] = Atom.x[0] - boxsize[0]

        # if Atom.x[1] < -0.5 * boxsize[1]: Atom.x[1] = Atom.x[1] + boxsize[1]
        # if Atom.x[1] >= 0.5 * boxsize[1]: Atom.x[1] = Atom.x[1] - boxsize[1]

        # if Atom.x[2] < -0.5 * boxsize[2]: Atom.x[2] = Atom.x[2] + boxsize[2]
        # if Atom.x[2] >= 0.5 * boxsize[2]: Atom.x[2] = Atom.x[2] - boxsize[2]

        Atom.x[0] = Atom.x[0] - np.floor(Atom.x[0] / boxsize[0]) * boxsize[0]
        Atom.x[1] = Atom.x[1] - np.floor(Atom.x[1] / boxsize[1]) * boxsize[1]
        Atom.x[2] = Atom.x[2] - np.floor(Atom.x[2] / boxsize[2]) * boxsize[2]

        # Atom.x[0] %= boxsize[0]
        # Atom.x[1] %= boxsize[1]
        # Atom.x[2] %= boxsize[2]

        Atom.v[0] += 0.5 * (Atom.a[0] + Atom.a_new[0]) * dt
        Atom.v[1] += 0.5 * (Atom.a[1] + Atom.a_new[1]) * dt
        Atom.v[2] += 0.5 * (Atom.a[2] + Atom.a_new[2]) * dt
        
        Atom.a[0] = Atom.a_new[0]
        Atom.a[1] = Atom.a_new[1]
        Atom.a[2] = Atom.a_new[2]

# Write coordinates to a .pdb trajectory file
def writeFrame(step, AtomList, boxsize):
    with open("traj.pdb", "a+") as file:
        
        file.write("MODEL        {}\n".format(step))
        file.write("CRYST1  {:7.3f}  {:7.3f}  {:7.3f}  90.00  90.00  90.00 P 1           1\n".format(boxsize[0], boxsize[1], boxsize[2]))
        
        for Atom in AtomList:
            file.write("{:6s}{:5d} {:^4s}{:1s}{:4s}{:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}\n".format('ATOM', Atom.idx, Atom.name, '', Atom.name, '', Atom.idx, '', Atom.x[0], Atom.x[1], Atom.x[2]))
        
        file.write('TER\nENDMDL\n')

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

    for step in range(nsteps):
        # 2 COMPUTE FORCE
        get_accelerations(AtomList, boxsize)

        # 3 UPDATE CONFIGURATION
        integrate(AtomList, dt, boxsize)

        # 4 OUTPUT STEP
        if (step % 100 == 0):
            writeFrame(step + 1, AtomList, boxsize)
        positionList.append([Atom.x[0] for Atom in AtomList])

    return positionList

# MAIN #########################################################################

AtomList = [
    Atom(1, 'ARG', [1.0, 0, 0], 39.948, 0, 0.0103, 3.4),
    Atom(2, 'ARG', [5.0, 0, 0], 39.948, 0, 0.0103, 3.4),
    Atom(3, 'ARG', [9.0, 0, 0], 39.948, 0, 0.0103, 3.4)
    ]

sim_pos = run_md(AtomList=AtomList, dt=0.1, nsteps=200000, T=300, boxsize=[30.0, 30.0, 30.0])

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
