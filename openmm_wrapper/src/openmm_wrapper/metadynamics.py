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

import logging
import openmm as mm
import openmm.app as app
import openmm_wrapper as my

import pathlib
import json
import copy
import numpy as np
import pprint as pp


def initialiseMetadynamics(setup, modeller, system):
    logger = logging.getLogger("dynamicEntropy")

    logger.critical("Metadynamcis ...")
    mtdCommands = setup.config["metadynamics"]
    dumpGaussianParameters(mtdCommands)

    collectiveVariables = []
    for key, cv in mtdCommands.items():
        if "cv" in key:
            if cv is None:
                continue

            logger.info("{}".format(key))
            dumpCVParameters(cv)
            collectiveVariablesForce = my.createCollectiveVariables(
                cv, modeller.topology, system
            )

            metadynamicsWalls(system, collectiveVariablesForce, cv)

            newCV = app.BiasVariable(
                collectiveVariablesForce,
                cv["lowerLimit"],
                cv["upperLimit"],
                cv["sigma"],
                periodic=cv["periodic"],
                gridWidth=cv["gridWidth"],
            )
            collectiveVariables.append(newCV)

    MTD = my.newMetadynamicsKernel(
        system,
        collectiveVariables,
        setup.config["md"]["temperature"],
        mtdCommands["biasFactor"],
        mtdCommands["height"],
        mtdCommands["frequency"],
        saveFrequency=mtdCommands["saveFrequency"],
        biasDir=mtdCommands["biasDir"],
    )

    return MTD


def metadynamicsWalls(system, cvForce, config):
    """
    Add a wall on the CV using an automatically defined function.
    It work only for nonbonded type of CVs.

    """
    logger = logging.getLogger("dynamicEntropy")

    wallString = config["expression"]
    if "uwall" in config:
        if config["uwall"] is not None:
            wallExpr = "step(cv)*k*cv^2; cv=({}{:+}); k={}".format(
                wallString, -config["uwall"], config["uwallK"]
            )
            x = wallExpr.split(";")
            logger.info("  {:40s} = {} ".format("Lower Wall expression", x[0]))
            for i in range(1, len(x)):
                logger.info("  {:40s} = {} ".format(" ", x[i].strip()))
            addWall(system, cvForce, wallExpr)

    if "lwall" in config:
        if config["lwall"] is not None:
            wallExpr = "(1-step(cv))*k*cv^2; cv=({}{:+}); k={}".format(
                wallString, -config["lwall"], config["lwallK"]
            )
            x = wallExpr.split(";")
            logger.info("  {:40s} = {} ".format("Lower Wall expression", x[0]))
            for i in range(1, len(x)):
                logger.info("  {:40s} = {} ".format(" ", x[i].strip()))
            addWall(system, cvForce, wallExpr)


def addWall(system, cvForce, wallExpr):
    wallForce = copy.copy(cvForce)
    wallForce.setEnergyFunction(wallExpr)
    wallForce.setName("MTDWall")

    # Set the exclustions as in the nonBondedForce
    if isinstance(wallForce, mm.CustomNonbondedForce):
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                for i in range(0, f.getNumExceptions()):
                    lst = f.getExceptionParameters(i)
                    wallForce.addExclusion(lst[0], lst[1])
    system.addForce(wallForce)


def runMetadynamics(setup, MTD, simulation):
    """
    Wrapper to run metadynamics steps and save the final free energy.
    """
    logger = logging.getLogger("dynamicEntropy")
    logger.critical("#--- Running metadynamics ------------------#")
    if setup.debug:
        my.computeEnergyFromContext(simulation.context)
        logger.info(
            "  {:40s} = {} ".format(
                "Collective Variables", MTD.getCollectiveVariables(simulation)
            )
        )
        return

    MTD.step(simulation, setup.config["md"]["numberOfSteps"])
    np.save("freeEnergy", MTD.getFreeEnergy())


def dumpCVParameters(config):
    logger = logging.getLogger("dynamicEntropy")
    for key, val in config.items():
        logger.info("  {:40s} = {} ".format(key, val))


def dumpGaussianParameters(config):
    logger = logging.getLogger("dynamicEntropy")
    for key, val in config.items():
        if "cv" in key:
            continue
        logger.info("  {:40s} = {} ".format(key, val))
