def createPlumed(plumed_file, system):
    from openmmplumed import PlumedForce

    with open(plumed_file, "r") as f:
        script = f.read()
    system.addForce(PlumedForce(script))
    return "PLUMED"
