### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import sys
import copy
import logging
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import new_openmm_wrapper as my

import pathlib


def initialiseMolecularDynamics(simulation, positions, temperature):
    my.pretty_log(title="Initialising molecular dynamics", sep=True)
    my.pretty_log(f" Setting temperature to {temperature}", indent=1)
    # Initialise positions and velocities
    # Need to make a function to deal with restarts
    simulation.context.setPositions(positions)
    simulation.context.setVelocitiesToTemperature(temperature, 12345)


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

    # temp hack
    # addException(particle1, particle2, chargeProd, sigma, epsilon, replace=False)
    for force in system.getForces():
        if isinstance(force, mm.NonbondedForce):
            for i in range(5):
                force.addException(i, 5, 0, 1, 0)
        if isinstance(force, mm.CustomNonbondedForce):
            for i in range(5):
                force.addExclusion(i, 5)

    simulation = my.createSimulation(
        modeller.topology,
        system,
        integrator,
        mm.Platform.getPlatformByName(setup.config["basic"]["Platform"]),
        setup.properties,
    )

    reporters = my.simulationReporters(
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
        simulation, modeller.positions, setup.config["md"]["temperature"]
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

    if config["restartOutput"]["file"] is not None:
        my.pretty_log(f"Writing restart file ({config['restartOutput']['file']})")
        # my.dumpRestartFiles(simulation, config["restartOutput"]["file"])


# to delete
def oldAddReporters(simulation, setup, modeller, system):
    """
    Adds all configured reporters to the simulation object.

    Parameters:
        simulation (Simulation): The simulation object to add reporters to.
        setup (Setup): The setup object containing configuration parameters.
        modeller (Modeller): The modeller object representing the system topology.
        system (System): The system object representing the molecular system.

    Returns:
        None
    """
    logger = logging.getLogger("dynamicEntropy")

    mdConfig = setup.config["md"]

    # Screen output reporter
    if mdConfig["screenOutput"] is not None:
        mdConfig["screenOutput"]["totalSteps"] = mdConfig["numberOfSteps"]
        if mdConfig["screenOutput"]["reportInterval"] == 0:
            mdConfig["screenOutput"]["reportInterval"] = int(
                mdConfig["numberOfSteps"] / 100
            )

        if mdConfig["simulationTime"] is not None:
            mdConfig["screenOutput"]["progress"] = False
            mdConfig["screenOutput"]["remainingTime"] = False

        if setup.debug:
            my.dumpInfo("screen output", mdConfig["screenOutput"])
        else:
            logger.info(
                "  {:40s} = {}".format(
                    "screen output frequency",
                    mdConfig["screenOutput"]["reportInterval"],
                )
            )

        # Find barostat if present and add to reporter
        for force in system.getForces():
            if isinstance(force, mm.MonteCarloBarostat):
                mdConfig["screenOutput"]["barostat"] = force
                mdConfig["logOutput"]["barostat"] = force
                break

        if mdConfig["screenOutput"]["file"] is not None:
            if mdConfig["screenOutput"]["file"].lower() == "stdout":
                mdConfig["screenOutput"]["file"] = sys.stdout
            elif mdConfig["screenOutput"]["file"].lower() == "stderr":
                mdConfig["screenOutput"]["file"] = sys.stderr

            simulation.reporters.append(
                my.ExtendedStateDataReporter(**mdConfig["screenOutput"])
            )

    # Log file reporter
    if mdConfig["logOutput"] is not None:
        if mdConfig["logOutput"]["reportInterval"] == 0:
            mdConfig["logOutput"]["reportInterval"] = mdConfig["reportInterval"]

        if setup.debug:
            my.dumpInfo("log output", mdConfig["logOutput"])
        else:
            logger.info(
                "  {:40s} = {}".format("log file", mdConfig["logOutput"]["file"])
            )
            logger.info(
                "  {:40s} = {}".format(
                    "log frequency", mdConfig["logOutput"]["reportInterval"]
                )
            )
        simulation.reporters.append(
            my.ExtendedStateDataReporter(**mdConfig["logOutput"])
        )

    # Trajectory output reporter
    if (
        mdConfig["trajectoryOutput"] is not None
        and mdConfig["trajectoryOutput"]["reportInterval"] >= 0
        and mdConfig["trajectoryOutput"]["file"] is not None
    ):
        if mdConfig["trajectoryOutput"]["reportInterval"] == 0:
            mdConfig["trajectoryOutput"]["reportInterval"] = mdConfig["reportInterval"]

        if setup.debug:
            my.dumpInfo("trajectory output", mdConfig["trajectoryOutput"])
        else:
            logger.info(
                "  {:40s} = {}".format(
                    "trajectory file", mdConfig["trajectoryOutput"]["file"]
                )
            )
            logger.info(
                "  {:40s} = {}".format(
                    "trajectory frequency",
                    mdConfig["trajectoryOutput"]["reportInterval"],
                )
            )

        file_extension = pathlib.Path(mdConfig["trajectoryOutput"]["file"]).suffix
        if file_extension == ".dcd":
            simulation.reporters.append(app.DCDReporter(**mdConfig["trajectoryOutput"]))
        elif file_extension == ".xtc":
            simulation.reporters.append(app.XTCReporter(**mdConfig["trajectoryOutput"]))
        elif file_extension == ".pdb":
            simulation.reporters.append(app.PDBReporter(**mdConfig["trajectoryOutput"]))
        else:
            raise Exception("Unknown trajectory format")

    # Osmotic pressure reporter
    if setup.config["osmotic"] is not None:
        value = setup.config["osmotic"]
        if all(isinstance(value, dict) for value in value.values()):
            for k1, v1 in value.items():
                simulation.reporters.append(
                    my.OsmoticPressureReporter(
                        setup.config["input"]["runID"],
                        k1,
                        v1,
                        mdConfig["temperature"],
                        modeller.topology,
                        simulation.context,
                    )
                )
        else:
            simulation.reporters.append(
                my.OsmoticPressureReporter(
                    setup.config["input"]["runID"],
                    "o",
                    value,
                    mdConfig["temperature"],
                    modeller.topology,
                    simulation.context,
                )
            )

    # Catalytic reaction reporter
    if setup.config["catalyst"] is not None:
        simulation.reporters.append(
            my.CatalyticReactionReporter(
                setup.config["catalyst"],
            )
        )

    # Force reporter for testing
    if setup.config["test_forces"]:
        simulation.reporters.append(my.forceReporter("forces.dat", 1000))

    # New FEP
    # if setup.config["newfep"]:
    simulation.reporters.append(my.newFepReporter(10, setup))

    # Custom COM remover reporter
    if mdConfig["CMMotionRemover"] is not None:
        if mdConfig["CMMotionRemover"].get("type") == "custom":
            settings = mdConfig["CMMotionRemover"]
            simulation.reporters.append(
                my.customCOMRemover(settings, modeller.topology)
            )
