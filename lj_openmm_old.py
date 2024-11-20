#!/usr/bin/env python
# coding: utf-8



import openmm as mm
from openmm import app
from openmm import unit
from openmm.app.pdbfile import PDBFile
import numpy as np
import sys


# Parameters
n_particles = 500  # Total number of particles
fraction_A = 0.5   # Fraction of particles of type A
box_size = 3.0
report_interval = 20
temperature = 300 * unit.kelvin
pressure = 10.0 * unit.atmosphere
simulation_steps = 100000  # Number of simulation steps
timestep = 2.0 * unit.femtoseconds

# Lennard-Jones parameters
sigma_A = 0.34 * unit.nanometer
epsilon_A = 0.5 * unit.kilojoule_per_mole
sigma_B = 0.47 * unit.nanometer
epsilon_B = 0.8 * unit.kilojoule_per_mole

# Mixing rules (Lorentz-Berthelot)
sigma_AB = (sigma_A + sigma_B) / 2
epsilon_AB = np.sqrt(epsilon_A * epsilon_B)

# Create particle types
n_A = int(fraction_A * n_particles)
n_B = n_particles - n_A
positions = np.random.uniform(low=-box_size/2, high=box_size/2, size=(n_particles, 3))

print("Shape of positions:", positions.shape)  # Should be (n_particles, 3)

# OpenMM System
system = mm.System()

# Add particles to the system
particle_types = ['A'] * n_A + ['B'] * n_B
for _ in range(n_particles):
    system.addParticle(1.0 * unit.dalton)  # Assign equal masses to all particles

print("Number of particles in system:", system.getNumParticles())
print("Number of positions provided:", len(positions))

# Define periodic box vectors
box_vectors = mm.Vec3(box_size * unit.nanometer, 0, 0), mm.Vec3(0, box_size * unit.nanometer, 0), mm.Vec3(0, 0, box_size * unit.nanometer)
system.setDefaultPeriodicBoxVectors(*box_vectors)

# Create a nonbonded force for Lennard-Jones interactions
nonbonded = mm.NonbondedForce()
nonbonded.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
nonbonded.setCutoffDistance(1.2 * sigma_B)
nonbonded.setUseDispersionCorrection(True)

# Assign Lennard-Jones parameters
for i, ptype in enumerate(particle_types):
    if ptype == 'A':
        nonbonded.addParticle(0.0, sigma_A, epsilon_A)  # q = 0.0 (neutral particles)
    else:
        nonbonded.addParticle(0.0, sigma_B, epsilon_B)
        
# Define cross-interactions explicitly for A-B pairs
for i in range(n_particles):
    for j in range(i + 1, n_particles):
        if particle_types[i] != particle_types[j]:  # A-B interaction
            nonbonded.addException(i, j, 0.0, sigma_AB, epsilon_AB)

system.addForce(nonbonded)

# Add a Monte Carlo barostat
# system.addForce(mm.MonteCarloBarostat(pressure, temperature))

integrator = mm.LangevinIntegrator(temperature, 1/unit.picosecond, timestep)

# Simulation setup
platform = mm.Platform.getPlatformByName('CUDA')  # Use CPU for simplicity

topology = app.Topology()
chain = topology.addChain()

# Define particle types and names
particle_types = ['A'] * n_A + ['B'] * n_B  # Define particle types
for i, ptype in enumerate(particle_types):
    residue = topology.addResidue("Particle", chain)
    atom_name = f"{ptype}_{i+1}"  # Example: "A_1", "B_2", etc.
    topology.addAtom(atom_name, app.Element.getByAtomicNumber(1), residue)




topology.setPeriodicBoxVectors(box_size * unit.nanometer * np.identity(3))
print (topology.getUnitCellDimensions())




simulation = app.Simulation(topology, system, integrator, platform)




# Set initial positions
simulation.context.setPositions(positions * unit.nanometer)




# Minimize energy
print("Minimizing energy...")
simulation.minimizeEnergy()

positions = simulation.context.getState(getPositions=True).getPositions()
with open("minimized.pdb", "w") as f:
    PDBFile.writeFile(simulation.topology, positions, f)

# Equilibrate
print("Equilibrating...")
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(10000)

# Add reporters
simulation.reporters.append(app.StateDataReporter(sys.stdout, report_interval, step=True,
                                                  potentialEnergy=True, temperature=True, density=True))
simulation.reporters.append(app.DCDReporter('trajectory.dcd', report_interval))  # Save trajectory
simulation.reporters.append(app.CheckpointReporter('checkpoint.chk', report_interval))  # Save checkpoints

# Production run
print("Running production...")
simulation.step(simulation_steps)
