
import numpy as np
import sys, itertools, yaml

import openmm as mm
from openmm import app
from openmm import unit
from openmm.app.pdbfile import PDBFile

class Sim(object):
    def __init__(self):
        """
        Class for running simulations of Lennard Jones mixtures in OpenMM

        """
        print ("\n# --------------------------------------------------")
        print ("# Generated instance of Sim class in lj_openmm.py")
        print ("# A Python program for the simulation of mixtures of")
        print ("# Lennard-Jones particles")
        print ("# --------------------------------------------------\n")

    def read_input(self, fileinp):
        """ Reads input from file 

        Parameters
        ----------
        fileinp : str
            Input file in yml format

        """
        try:
            with open(fileinp) as stream:
                try:
                    config = yaml.safe_load(stream)
                except yaml.YAMLError as e:
                    print (e)
        except FileNotFoundError as e:
            print ("\nFile %s not found"%fileinp)
            print (e)   
            sys.exit()

        # Get config parameters for the simulation
        for k,v in config.items():
            setattr(self, k, v)

    def make_system(self):
        """
        Add particles to the system

        """
        system = mm.System()

        # Add particles to the system
        n_particles = 0
        for k,v in self.components.items():
            n_particles += v['nmol'] 
            for _ in range(v['nmol']):
                system.addParticle(v['mass'] * unit.dalton)
        print("# Number of particles in system:", system.getNumParticles())

        # Define periodic box vectors
        box_size = self.box
        box_vectors = mm.Vec3(box_size[0] * unit.nanometer, 0, 0),\
                    mm.Vec3(0, box_size[1] * unit.nanometer, 0), \
                    mm.Vec3(0, 0, box_size[2] * unit.nanometer)
        system.setDefaultPeriodicBoxVectors(*box_vectors)
        self.system = system

    def make_topology(self):
        """
        Creates topology

        """
        topology = app.Topology()
        chain = topology.addChain()

        # Define particle types and names
        for k,v in self.components.items():
            residue = topology.addResidue("Particle", chain)
            for i in range(v['nmol']):
                atom_name = f"{k}{i+1}"  # Example: "A_1", "B_2", etc.
                topology.addAtom(atom_name, app.Element.getByAtomicNumber(1), residue)

        # Generate random positions
        box_size = self.box
        topology.setPeriodicBoxVectors(box_size * unit.nanometer * np.eye(3))

        self.topology = topology

    def force_field(self):
        """
        Generates force field

        """
        # Lennard-Jones parameters
        for k, v in self.components.items():
            self.components[k]['sigma'] *= unit.nanometer
            self.components[k]['epsilon'] *= unit.kilojoule_per_mole 

        # Create a nonbonded force for Lennard-Jones interactions
        nonbonded = mm.NonbondedForce()
        nonbonded.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
        nonbonded.setCutoffDistance(1.5 * unit.nanometer)
        nonbonded.setUseDispersionCorrection(True)

        # Assign Lennard-Jones parameters
        for k,v in self.components.items():
            for i in range(v['nmol']):
                nonbonded.addParticle(0.0, v['sigma'] , v['epsilon'])

        # Define cross-interactions explicitly for A-B pairs
        if len(self.components) > 1:
            print ("# Defining cross interaction terms")
            allkeys = self.components.keys()
            for k,kk in itertools.combinations(allkeys, 2):
                print (k,kk)
                # Mixing rules (Lorentz-Berthelot)
                sigma = (self.components[k]['sigma'] + self.components[kk]['sigma'])/2
                epsilon = np.sqrt(self.components[k]['epsilon'] * self.components[kk]['epsilon'])
                print (sigma, epsilon)
                for i in range(self.components[k]['nmol']):
                    for j in range(self.components[kk]['nmol']):
                        jj = self.components[k]['nmol'] + j
                        try:
                            nonbonded.addException(i, jj, 0.0, sigma, epsilon)
                        except Exception as e:
                            print (e)
                            sys.exit()
        self.system.addForce(nonbonded)
        print ("# Generated force field")

    def make_simulation(self):
        """
    
        """
        # Simulation setup
        try:
            platform = mm.Platform.getPlatformByName('CUDA')
            print ("# Running on GPU")
        except Exception as e:
            platform = mm.Platform.getPlatformByName('CPU')
            print ("# CUDA Platform not available; running on CPU")

        if self.pressure:
            self.system.addForce(mm.MonteCarloBarostat(self.pressure * unit.bar,\
                                self.temperature * unit.kelvin))
            
        integrator = mm.LangevinIntegrator(self.temperature * unit.kelvin, \
                                1/unit.picosecond, self.timestep)
        simulation = app.Simulation(self.topology, self.system, integrator, platform)

        # Set initial positions
        box_size = self.box
        min_box_size = np.min(box_size)
        n_particles = sum([v['nmol'] for k,v in self.components.items()])
        positions = np.random.uniform(low=-min_box_size/2, high=min_box_size/2, \
                                          size=(n_particles, 3))
        print("# Number of positions provided:", len(positions))
        simulation.context.setPositions(positions * unit.nanometer)
        self.simulation = simulation

    def run(self):
        """

        """
        # Minimize energy
        print("# Minimizing energy...")
        self.simulation.minimizeEnergy()

        positions = self.simulation.context.getState(getPositions=True).getPositions()
        froot = self.output_name
        with open(froot + "_min.pdb", "w") as f:
            PDBFile.writeFile(self.simulation.topology, positions, f)

        self.simulation.reporters.append(app.StateDataReporter(froot + "_equil.csv", \
                        100, time=True, potentialEnergy=True, temperature=True, \
                        density=True))
        self.simulation.reporters.append(app.DCDReporter(froot + \
                            '_equil.dcd', 100))

        if self.pressure:
            # Equilibration in NPT
            self.simulation.context.setVelocitiesToTemperature(self.temperature)
            self.simulation.context.setPositions(positions)
            print("# Equilibrating...")
            self.simulation.step(self.equilibration_steps)

            positions = self.simulation.context.getState(getPositions=True).getPositions()
            with open(froot + "_equil.pdb", "w") as f:
                PDBFile.writeFile(self.simulation.topology, positions, f)

        # Remove the MC barostat
        if self.pressure:
            forces = self.simulation.system.getForces()
            for c, f in enumerate(forces):
                if isinstance(f, mm.MonteCarloBarostat):
                    self.simulation.system.removeForce(c)


        # Production run in NVT
        print("# Production run...")

        self.simulation.reporters = []
        self.simulation.reporters.append(app.StateDataReporter(froot + "_prod.csv", 
                            self.report_interval, time=True, potentialEnergy=True, \
                            temperature=True, density=True))
        self.simulation.reporters.append(app.DCDReporter(froot + '_prod.dcd', \
                            self.report_interval))

        self.simulation.step(self.production_steps)
        positions = self.simulation.context.getState(getPositions=True).getPositions()
        with open(froot + "_final.pdb", "w") as f:
            PDBFile.writeFile(self.simulation.topology, positions, f)