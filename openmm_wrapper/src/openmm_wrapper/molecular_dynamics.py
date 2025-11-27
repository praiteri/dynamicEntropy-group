### Copyright (C) 2023  Paolo Raiteri

### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
### GNU General Public License for more details.

### You should have received a copy of the GNU General Public License
### along with this program.  If not, see <https://www.gnu.org/licenses/>.

import sys
import copy
import logging
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import openmm_wrapper as my

import numpy as np
import pathlib


def initialiseMolecularDynamics(setup, modeller, system):
    """
    Initializes the molecular dynamics simulation using the provided setup, modeller, and system.

    Parameters:
        setup (Setup): The setup object containing the configuration parameters for the simulation.
        modeller (Modeller): The modeller object representing the system topology.
        system (System): The system object representing the molecular system.

    Returns:
        Simulation: The initialized simulation object.

    Raises:
        None

    """
    logger = logging.getLogger("dynamicEntropy")

    if system.getNumParticles() == 1:
        logger.warning("----------------------------------#")
        logger.warning(" only one atom deleting CMMotionRemover")
        logger.warning("----------------------------------#")
        for i, force in enumerate(system.getForces()):
            if isinstance(force, mm.CMMotionRemover):
                system.removeForce(i)

    # Screen output
    setup.dumpParametersMD()

    # Configuration parameters for MD
    mdConfig = copy.deepcopy(setup.config["md"])
    mdConfig.pop("simulationTime")

    # Create integrator object
    integrator = my.integrator(system, mdConfig)

    # Create simulatin object
    simulation = my.createSimulation(
        modeller.topology,
        system,
        integrator,
        mm.Platform.getPlatformByName(setup.config["basic"]["Platform"]),
        setup.properties,
    )

    # simulation.reporters.append(my.removeTorque(system, 1000))

    # restart old simulation
    if mdConfig["restartFrom"] is not None:
        logger.critical("Restarting MD simulation ...")
        logger.info("  {:40s} = {}".format("file", mdConfig["restartFrom"]))

        simulation.context.setStepCount(0)

        file_extension = pathlib.Path(mdConfig["restartFrom"]).suffix
        if file_extension == ".xml":
            simulation.loadState(mdConfig["restartFrom"])
            simulation.context.reinitialize(preserveState=True)

        elif file_extension == ".chk":
            simulation.loadCheckpoint(mdConfig["restartFrom"])

        elif file_extension == ".pdb":
            atomic_coordinates, cell = my.read_pdb_coordinates(
                mdConfig["restartFrom"], unit="nm"
            )

            pos = np.array([a["coord"] for a in atomic_coordinates])

            hmat = app.internal.unitcell.computePeriodicBoxVectors(*cell)

            simulation.context.setPositions(pos)
            simulation.context.setPeriodicBoxVectors(*hmat)

            simulation.context.setVelocitiesToTemperature(
                mdConfig["temperature"],
                setup.rng.integers(low=1, high=99999, size=1)[0],
            )

            # my.computeEnergyFromContext(simulation.context)
            # quit()

    else:
        # Initialise positions and velocities
        simulation.context.setPositions(modeller.positions)

        simulation.context.setVelocitiesToTemperature(
            mdConfig["temperature"], setup.rng.integers(low=1, high=99999, size=1)[0]
        )

    if setup.config["md"]["minimise"]:
        logger.critical("Energy minimisation (MD) ...")
        # reporter = my.minimizationReporter()
        simulation.minimizeEnergy()
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        pos = state.getPositions(asNumpy=True)
        simulation.context.reinitialize()
        simulation.context.setPositions(pos)
        simulation.context.setVelocitiesToTemperature(
            setup.config["md"]["temperature"],
            setup.rng.integers(low=1, high=99999, size=1)[0],
        )

    # Screen output
    if mdConfig["screenOutput"] is not None:
        mdConfig["screenOutput"]["totalSteps"] = setup.config["md"]["numberOfSteps"]
        if mdConfig["screenOutput"]["reportInterval"] == 0:
            mdConfig["screenOutput"]["reportInterval"] = int(
                setup.config["md"]["numberOfSteps"] / 100
            )

        if setup.config["md"]["simulationTime"] is not None:
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

        if mdConfig["screenOutput"]["file"] is not None:
            if mdConfig["screenOutput"]["file"].lower() == "stdout":
                mdConfig["screenOutput"]["file"] = sys.stdout
            elif mdConfig["screenOutput"]["file"].lower() == "stderr":
                mdConfig["screenOutput"]["file"] = sys.stderr
            simulation.reporters.append(
                app.StateDataReporter(**mdConfig["screenOutput"])
            )

    # Log file
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
        simulation.reporters.append(app.StateDataReporter(**mdConfig["logOutput"]))

    # Trajectory output
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

    # Osmotic pressure
    if setup.config["osmotic"] is not None:
        value = setup.config["osmotic"]
        if all(isinstance(value, dict) for value in value.values()):
            for k1, v1 in value.items():
                simulation.reporters.append(
                    my.OsmoticPressureReporter(
                        setup.config["input"]["runID"],
                        k1,
                        v1,
                        setup.config["md"]["temperature"],
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
                    setup.config["md"]["temperature"],
                    modeller.topology,
                    simulation.context,
                )
            )

    # Catalytic reaction
    if setup.config["catalyst"] is not None:
        simulation.reporters.append(
            my.CatalyticReactionReporter(
                setup.config["catalyst"],
            )
        )

        simulation.reporters.append(
            my.CatalyticReactionReporter("forces.dat", simulation.context)
        )

    # simulation.reporters.append(
    #     my.forceReporter("forces.dat",1000)
    # )

    return simulation


def runMolecularDynamics(setup, simulation):
    """
    Runs molecular dynamics simulation.

    Args:
        setup (Setup): The setup object containing configuration settings.
        simulation (Simulation): The simulation object.

    Returns:
        None
    """
    if setup.debug:
        my.computeEnergyFromContext(simulation.context)

    logger = logging.getLogger("dynamicEntropy")
    logger.critical("#--- Running molecular dynamics ------------#")

    if setup.config["md"]["simulationTime"] is not None:
        h, m, s = [0, 0, 0]
        if isinstance(setup.config["md"]["simulationTime"], int):
            t = setup.config["md"]["simulationTime"]
        else:
            string = setup.config["md"]["simulationTime"]
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

        simulation.runForClockTime(t)
    else:
        simulation.step(int(setup.config["md"]["numberOfSteps"]))

    if setup.config["md"]["restartOutput"]["file"] is not None:
        logger.critical("#--- Writing restart file ------------------#")
        my.dumpRestartFiles(simulation, setup.config["md"]["restartOutput"]["file"])
