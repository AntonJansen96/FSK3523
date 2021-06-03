#!/bin/python3

import numpy as np
import matplotlib.pyplot as plt

class constants:
    R   = 8.3144621       # J/mol/K
    N   = 6.02214129e23   # 1/mol
    k_b = 1.380649e-23    # J/K
    e   = 1.602176634e-19 # C
    u   = 1.660539066e-27 # kg

class Atom:
    def __init__(self, idx, name, x, mass, charge, lj_epsilon, lj_sigma):
        self.idx        = idx
        self.name       = name
        self.mass       = mass
        self.charge     = charge
        self.lj_epsilon = lj_epsilon
        self.lj_sigma   = lj_sigma
        self.x          = x
        self.v          = 0
        self.a          = 0
        self.a_new      = 0

def LJ_force(r, epsilon, sigma):
    repulsive  = 48 * epsilon * sigma**12 / r**13
    attractive = 24 * epsilon * sigma**6  / r**7

    return repulsive - attractive

def generate_velocities(AtomList, T):
    for Atom in AtomList:
        BoltzmannFactor = constants.k_b * T
        BoltzmannFactor = BoltzmannFactor / (Atom.mass * constants.e)
        BoltzmannFactor = BoltzmannFactor**0.5
        Atom.v          = BoltzmannFactor * (np.random.rand() - 0.5)

def get_accelerations(AtomList):
    def reduce(forcegrid):
        reduced = len(forcegrid) * [0]

        for i in range(0, len(forcegrid)):
            for j in range(0, len(forcegrid)):
                reduced[i] += forcegrid[j][i]

        return reduced
    
    accel_x = [[0] * len(AtomList) for i in range(len(AtomList))]
        
    for i in range(0, len(AtomList) - 1):
        for j in range(i + 1, len(AtomList)):
            
            r_x = AtomList[j].x - AtomList[i].x
            
            rmag = (r_x**2)**0.5
            
            force_scalar = LJ_force(rmag, AtomList[i].lj_epsilon, AtomList[i].lj_sigma)

            force_x = force_scalar * r_x / rmag
            
            accel_x[i][j] =   force_x / AtomList[i].mass
            accel_x[j][i] = - force_x / AtomList[j].mass

    accels = reduce(accel_x)
    
    for i in range(0, len(AtomList)):
        AtomList[i].a_new = accels[i]

# Velocity-Verlet integrator.
def integrate(AtomList, dt):
    for Atom in AtomList:
        Atom.x += Atom.v * dt + 0.5 * Atom.a * dt * dt
        Atom.v += 0.5 * (Atom.a + Atom.a_new) * dt
        Atom.a  = Atom.a_new

# Write coordinates to a .pdb trajectory file
def writeFrame(AtomList, step):
    with open("traj.pdb", "a+") as file:
        
        file.write("MODEL        {}\n".format(step))
        
        for Atom in AtomList:
            file.write("{:6s}{:5d} {:^4s}{:1s}{:4s}{:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}\n".format('ATOM', Atom.idx, Atom.name, '', Atom.name, '', Atom.idx, '', Atom.x, 0.0, 0.0))
        
        file.write('TER\nENDMDL\n')

# Run MD simulation.
def run_md(AtomList, dt, nsteps, T):
    # 1 INPUT INITIAL CONDITIONS

        # Initialize list of positions for old output function.
    positionList = [[Atom.x for Atom in AtomList]]
    
        # Initial velocity from Maxwell-Boltzmann distribution.
    generate_velocities(AtomList, T)

        # Update user.
    print("initial positions:  {:.4f}, {:.4f}, {:.4f}".format(AtomList[0].x, AtomList[1].x, AtomList[2].x))
    print("initial velocities: {:.4f}, {:.4f}, {:.4f}".format(AtomList[0].v, AtomList[1].v, AtomList[2].v))

    for step in range(nsteps):
        # 2 COMPUTE FORCE
        get_accelerations(AtomList)

        # 3 UPDATE CONFIGURATION
        integrate(AtomList, dt)

        # 4 OUTPUT STEP
        writeFrame(AtomList, step + 1)
        positionList.append([Atom.x for Atom in AtomList])

    return positionList

# MAIN #########################################################################

AtomList = [
    Atom(1, 'ARG', 1.0,  39.948, 0, 0.0103, 3.4),
    Atom(2, 'ARG', 5.0,  39.948, 0, 0.0103, 3.4),
    Atom(3, 'ARG', 10.0, 39.948, 0, 0.0103, 3.4)
    ]

sim_pos = run_md(AtomList=AtomList, dt=0.1, nsteps=10000, T=300)

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
