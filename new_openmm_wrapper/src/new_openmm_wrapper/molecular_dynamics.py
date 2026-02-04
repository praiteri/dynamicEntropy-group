### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import os
import sys
import openmm as mm
import openmm.unit as unit
import new_openmm_wrapper as my


def initialiseMolecularDynamics(simulation, positions, config, seed=12345):
    my.pretty_log(title="Initialising molecular dynamics", sep=True)
    simulation.context.setPositions(positions)

    temperature = config["temperature"]
    if config["restartFrom"] is None:
        # Initialise positions and velocities - no restartFrom
        my.pretty_log(f" Creating velocities for temperature = {temperature}", indent=1)
        simulation.context.setVelocitiesToTemperature(temperature, seed)

    else:
        if os.path.splitext(config["restartFrom"])[-1] == ".pdb":
            my.pretty_log(
                f" Restarting simulation from positions only {config["restartFrom"]}",
                indent=1,
            )
            my.pretty_log(
                f" Creating velocities for temperature = {temperature}", indent=1
            )
            simulation.context.setVelocitiesToTemperature(temperature, seed)

        elif os.path.splitext(config["restartFrom"])[-1] == ".xml":
            my.pretty_log(
                f" Restarting simulation from state file {config["restartFrom"]}",
                indent=1,
            )
            simulation.loadState(config["restartFrom"])

        elif os.path.splitext(config["restartFrom"])[-1] == ".chk":
            my.pretty_log(
                f" Restarting simulation from checkpoint {config["restartFrom"]}",
                indent=1,
            )
            simulation.loadCheckpoint(config["restartFrom"])

        else:
            my.pretty_log(
                f"Unknown restart file type, us (xml, chk or pdb)", logger="error"
            )
            sys.exit(1)


def molecular_dynamics(config):
    setup = my.simulationSetup(config)
    modeller, system = my.createSystem(setup)

    # creating integrator
    integrator = my.createIntegrator(system, setup.config["md"])

    if config.get("plumed", None) is not None:
        my.createPlumed(config["plumed"], system)

    # Create metadynamics object if requested
    # CV force needs to be created before the simulation object
    if config.get("metadynamics", None) is not None:
        my.pretty_log(title="Creating metadynamics object", sep=True)
        MTD = my.initialiseMetadynamics(setup, modeller, system)

    simulation = my.createSimulation(
        modeller.topology,
        system,
        integrator,
        mm.Platform.getPlatformByName(setup.config["basic"]["Platform"]),
        setup.properties,
    )

    _ = my.simulationReporters(
        simulation,
        runID=setup.config["input"]["runID"],
        reportInterval=setup.config["md"]["reportInterval"],
        numberOfSteps=setup.config["md"].get("numberOfSteps", None),
        configReporters=setup.config.get("reporters", None),
    )

    # Add custom reporters here
    # - COM remover
    # - Osmotic pressure
    # - Catalysis
    # - test forces
    # - ...

    initialiseMolecularDynamics(
        simulation,
        modeller.positions,
        setup.config["md"],
        seed=setup.rng.integers(low=1, high=99999, size=1)[0],
    )

    if config.get("metadynamics", None) is not None:
        my.pretty_log(title="Adding HILLS reporter")
        simulation.reporters.append(
            my.HILLSReporter("HILLS", setup.config["metadynamics"]["frequency"], MTD)
        )
        my.runMetadynamics(setup, MTD, simulation)

    else:
        runMolecularDynamics(setup.config["md"], simulation)


def runMolecularDynamics(config, simulation):
    """
    Runs molecular dynamics simulation.
    """

    if config["simulationTime"] is not None:
        h, m, s = [0, 0, 0]
        if isinstance(config["simulationTime"], int):
            t = config["simulationTime"] * unit.seconds
        else:
            string = config["simulationTime"]
            nfields = string.count(":")
            if nfields == 2:
                h, m, s = string.split(":")
            elif nfields == 1:
                h, m = string.split(":")
            elif nfields == 0:
                h = string.split(":")[0]
            else:
                raise Exception("Cannot parse simulation time")

            h = float(h) * unit.hours
            m = float(m) * unit.minutes
            s = float(s) * unit.seconds
            t = h + m + s

        my.pretty_log(title=f"Running molecular dynamics for {t} seconds", sep=True)
        simulation.runForClockTime(t)
    else:
        my.pretty_log(
            title=f"Running molecular dynamics for {config['numberOfSteps']} steps",
            sep=True,
        )
        simulation.step(int(config["numberOfSteps"]))
