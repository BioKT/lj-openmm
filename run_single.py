from openmm import unit
import lj_openmm

simulation = lj_openmm.Sim()

simulation.read_input("single.yml")

simulation.make_system()

simulation.make_topology()

simulation.force_field()

simulation.make_simulation()

simulation.run()
