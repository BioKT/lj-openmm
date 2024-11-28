from openmm import unit
import lj_openmm

for i in ["argon", "krypton", "xenon"]:
    simulation = lj_openmm.Sim()

    simulation.read_input("data/%s_slab.yml"%i)

    simulation.make_system()

    simulation.make_topology()

    simulation.force_field()

    simulation.make_simulation()

    simulation.run()
