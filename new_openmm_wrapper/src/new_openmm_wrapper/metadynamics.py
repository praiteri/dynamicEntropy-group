### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import logging
import copy
import numpy as np

import openmm as mm
import openmm.app as app
import new_openmm_wrapper as my


logger = logging.getLogger("dynamicEntropy")


def initialiseMetadynamics(setup, modeller, system):
    logger.critical("Metadynamcis ...")
    mtdCommands = setup.config["metadynamics"]

    log_dict = copy.copy(mtdCommands)
    log_dict.pop("cv0", None)
    log_dict.pop("cv1", None)
    log_dict.pop("cv2", None)
    my.pretty_log(title="Metadynamics parameters", sep=True, data=log_dict)

    collectiveVariables = []
    for key, cv in mtdCommands.items():
        if "cv" in key:
            if cv is None:
                continue

            collectiveVariablesForce = my.createCollectiveVariables(cv, modeller.topology, system)
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

    MTD = app.Metadynamics(
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

    wallString = config["expression"]
    if "uwall" in config:
        if config["uwall"] is not None:
            wallExpr = "step(cv)*k*cv^2; cv=({}{:+}); k={}".format(
                wallString, -config["uwall"], config["uwallK"]
            )
            x = wallExpr.split(";")
            addWall(system, cvForce, wallExpr)
            my.pretty_log(x, title="Defining upper wall:", indent=2)

    if "lwall" in config:
        if config["lwall"] is not None:
            wallExpr = "(1-step(cv))*k*cv^2; cv=({}{:+}); k={}".format(
                wallString, -config["lwall"], config["lwallK"]
            )
            x = wallExpr.split(";")
            logger.info("  {:40s} = {} ".format("Lower Wall expression", x[0]))
            addWall(system, cvForce, wallExpr)
            my.pretty_log(x, title="Defining lowpper wall:", indent=2)


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
    my.pretty_log(
        title=f"Running metadynamics for {setup.config['md']['numberOfSteps']} steps", sep=True
    )

    MTD.step(simulation, setup.config["md"]["numberOfSteps"])
    np.save("freeEnergy", MTD.getFreeEnergy())
    if setup.config["md"]["restartOutput"]["file"] is not None:
        my.pretty_log(f"Writing restart file ({setup.config['md']['restartOutput']['file']})")
        # my.dumpRestartFiles(simulation, setup.config["md"]["restartOutput"]["file"])


def createPlumed(plumed_file, system):
    from openmmplumed import PlumedForce

    with open(plumed_file, "r") as f:
        script = f.read()
    system.addForce(PlumedForce(script))
    return
