import sys
import yaml
import lj_openmm
#test
for system in ["neon", "methane", "krypton", "xenon"]:
    template = "data/%s_L20.yml"%system
    for t in range(10,50,10):
        with open(template, "r") as f:
            try:
                inp = yaml.safe_load(f)
                inp['temperature'] = t
                inp['output_name'] = "data/%s_n50_L20_T%g"%(system,t)
            except yaml.YAMLError as e:
                print (e)
                sys.exit()
        with open("data/%s_n50_L20_new.yml"%system, "w") as f:
            yaml.dump(inp, f)

        simulation = lj_openmm.Sim()
        simulation.read_input("data/%s_n50_L20_new.yml"%system)
        simulation.make_system()
        simulation.make_topology()
        simulation.force_field()
        simulation.make_simulation()
        simulation.run()
