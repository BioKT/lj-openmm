from openmm import unit
import lj_openmm

simulation = lj_openmm.Sim()

simulation.read_input("single.yml")

simulation.generate_system()

#simulation.force_field()