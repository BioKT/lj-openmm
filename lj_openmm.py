import sys, yaml

import openmm as mm
from openmm import app
from openmm import unit
from openmm.app.pdbfile import PDBFile
import numpy as np

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

    def generate_system(self):
        """
        Add particles to the system

        """
        system = mm.System()

        # Add particles to the system
        n_particles = 0
        for k,v in self.components.items():
            n_particles += v['nmol'] 
            for _ in range(v['nmol']):
                system.addParticle(1.0 * unit.dalton)  # Assign equal masses to all particles

        # Generate random positions
        box_size = self.box
        min_box_size = np.min(box_size)
        positions = np.random.uniform(low=-min_box_size/2, high=min_box_size/2, \
                                      size=(n_particles, 3))
        print("Number of particles in system:", system.getNumParticles())
        print("Number of positions provided:", len(positions))

        # Define periodic box vectors
        box_vectors = mm.Vec3(box_size[0] * unit.nanometer, 0, 0),\
                    mm.Vec3(0, box_size[1] * unit.nanometer, 0), \
                    mm.Vec3(0, 0, box_size[2] * unit.nanometer)
        system.setDefaultPeriodicBoxVectors(*box_vectors)

        self.system = system

    def force_field(self):
        """
        Generates force field

        """
        n_components = len(self.components)

        # Lennard-Jones parameters
        for k, v in self.components:
            self.components[k]['sigma'] *= unit.nanometer
            self.components[k]['epsilon'] *= unit.kilojoule_per_mole

#        # Mixing rules (Lorentz-Berthelot)
#        sigma_AB = (sigma_A + sigma_B) / 2
#        epsilon_AB = np.sqrt(epsilon_A * epsilon_B)        
#
#        # Create a nonbonded force for Lennard-Jones interactions
        nonbonded = mm.NonbondedForce()
        nonbonded.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
        nonbonded.setCutoffDistance(1.2 * unit.kilojoule_per_mole)
        nonbonded.setUseDispersionCorrection(True)

        # Assign Lennard-Jones parameters
        for k, v in self.components:
            nonbonded.addParticle(0.0, v['sigma'] , v['epsilon'])  

#        # Define cross-interactions explicitly for A-B pairs
#        for i in range(n_components):
#            for j in range(i + 1, n_components):
#                if particle_types[i] != particle_types[j]:  # A-B interaction
#                    nonbonded.addException(i, j, 0.0, sigma_AB, epsilon_AB)
#
        self.system.addForce(nonbonded)