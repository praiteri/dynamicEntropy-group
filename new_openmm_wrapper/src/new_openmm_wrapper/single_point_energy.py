### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import new_openmm_wrapper as my
import openmm as mm


def single_point_energy(config):
    """Single point energy"""

    setup = my.simulationSetup(config)
    modeller, system = my.createSystem(setup)

    simulation = my.createSimulation(
        modeller.topology,
        system,
        mm.VerletIntegrator(0),
        mm.Platform.getPlatformByName(setup.config["basic"]["Platform"]),
        setup.properties,
    )
    simulation.context.setPositions(modeller.positions)

    # my.pretty_log(title="Single point energy calculation", sep=True)
    my.computeEnergyFromContext(simulation.context)
