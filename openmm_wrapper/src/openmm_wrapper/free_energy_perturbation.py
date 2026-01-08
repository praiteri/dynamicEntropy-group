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
import openmm.unit as unit
import openmm_wrapper as my

from .free_energy_perturbation_reporter import *


def free_energy_perturbation(setup, modeller, system):
    """
    Perform free energy perturbation calculation using OpenMM.

    Parameters:
        - setup (object): The setup object containing configuration settings.
        - modeller (object): The modeller object.
        - system (object): The system object.

    Returns:
        None
    """
    logger = logging.getLogger("dynamicEntropy")
    logger.critical("#--- Free energy perturbation -------------#")

    # Create a list of systems for convenience
    _, s1 = my.createSystem(setup, first=False)
    _, s2 = my.createSystem(setup, first=False)
    listOfSystems = [system, s1, s2]

    # Initialise FEP variables
    lambdaValues, fepAtoms = initialise_FEP(setup, modeller, system)

    if setup.config["fep"]["type"].lower() in ["vdw", "ff"]:
        for s in listOfSystems:
            if setup.config["forcefield"]["isAMOEBA"]:
                vdwAlchemyForAMOEBA(s, fepAtoms)
            else:
                vdwAlchemy(s, fepAtoms)

    elif setup.config["fep"]["type"].lower() not in ["coul", "ff", "fep"]:
        raise Exception("Unknown FEP type: {}".format(setup.config["fep"]["type"]))

    # Scale the charge on the perturbed atoms
    for il, s in enumerate(listOfSystems):
        my.scaleCharges(s, fepAtoms, scale=lambdaValues["coul"][il])

    # Create the simulation object for the MD
    simulation = my.initialiseMolecularDynamics(setup, modeller, system)

    # Create a list of contexts fro FEP calculation
    listOfContexts = [
        simulation.context,
        my.createBasicContext(listOfSystems[1], setup),
        my.createBasicContext(listOfSystems[2], setup),
    ]

    # Scale the global parameters in the context
    for i, c in enumerate(listOfContexts):
        if setup.config["forcefield"]["isAMOEBA"]:
            c.setParameter("AmoebaVdwLambda", lambdaValues["vdw"][i])
        else:
            c.setParameter("lambdaCoul", lambdaValues["coul"][i])
            c.setParameter("lambdaFF", lambdaValues["ff"][i])
            c.setParameter("lambdaVdw", lambdaValues["vdw"][i])

        try:
            c.setParameter("lambdaFEP", lambdaValues["fep"][i])
        except:
            pass

    # If debugging just print out the energies
    if setup.config["fep"]["test"]:
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        pos = state.getPositions()
        box = state.getPeriodicBoxVectors()
        for c in listOfContexts:
            c.setPositions(pos)
            c.setPeriodicBoxVectors(box[0], box[1], box[2])
            my.computeEnergyFromContext(c)
        return

    # Energy minimisation
    # if setup.config["fep"]["minimise"]:
    #     logger.critical("Energy minimisation (FEP) ...")
    #     reporter = my.minimizationReporter()
    #     simulation.minimizeEnergy(reporter=reporter)
    #     state = simulation.context.getState(getPositions=True,getEnergy=True)
    #     pos = state.getPositions(asNumpy=True)
    #     # simulation.context.reinitialize()
    #     simulation.context.setPositions(pos)
    #     simulation.context.setVelocitiesToTemperature(
    #         setup.config["md"]["temperature"] , setup.rng.integers(low=1,high=99999,size=1)[0] )

    # Append FEP reporter
    # only the contexts for the perturbed states are required
    simulation.reporters.append(
        my.FEPReporter(setup, listOfContexts[1:], cc=listOfContexts[0])
    )

    # Run simulation
    logger.critical("#--- Running molecular dynamics ------------#")
    simulation.step(int(setup.config["md"]["numberOfSteps"]))

    if setup.config["md"]["restartOutput"]["file"] is not None:
        logger.critical("#--- Writing restart file ------------------#")
        my.dumpRestartFiles(
            simulation, setup.config["md"]["restartOutput"]["file"], XML=False
        )

    return


def initialise_FEP(setup, modeller, system):
    """
    Initializes the Free Energy Perturbation (FEP) calculation.

    Args:
        setup (object): The setup object containing configuration settings.
        modeller (object): The modeller object.
        system (object): The system object.

    Returns:
        tuple: A tuple containing the lambda values and the atoms involved in the FEP calculation.
    """
    logger = logging.getLogger("dynamicEntropy")
    fepType = setup.config["fep"]["type"].lower()
    # if fepType not in ["coul", "vdw"]:
    #     fepType = "ff"
    logger.info("  {:40s} = {} ".format("FEP type", fepType))

    # Get the lambda values
    lambdaValues = {
        "vdw": [1.0, 1.0, 1.0],
        "coul": [1.0, 1.0, 1.0],
        "ff": [0.0, 0.0, 0.0],
        "fep": [0.0, 0.0, 0.0],
    }

    if isinstance(setup.config["fep"]["lambda"], str):
        setup.config["fep"]["lambda"] = [
            float(x) for x in setup.config["fep"]["lambda"].split(",")
        ]
    lambdaValues[fepType] = setup.config["fep"]["lambda"]

    fepAtoms = []
    if fepType != "fep":
        if fepType == "vdw":
            lambdaValues["coul"] = [0.0, 0.0, 0.0]

        # Identify the atoms involved in the FEP
        fepAtoms = my.getAtoms(setup.config["fep"], modeller.topology)
        assert len(fepAtoms) > 0, "No atoms selected for FEP"
        logger.info("  {:40s} = {} ".format("Nr of atoms selected", len(fepAtoms)))
        if setup.debug:
            logger.info(" FEP: List of selected atoms")
            for i in range(0, len(fepAtoms), 10):
                logger.debug("FEP: {}".format(fepAtoms[i : i + 10]))

    logger.info(
        "  {:40s} = {} ".format("lambda values (coulomb)", lambdaValues["coul"])
    )
    logger.info("  {:40s} = {} ".format("lambda values (vdw)", lambdaValues["vdw"]))
    logger.info("  {:40s} = {} ".format("lambda values (ff)", lambdaValues["ff"]))
    logger.info("  {:40s} = {} ".format("lambda values (fep)", lambdaValues["fep"]))

    logger.info(
        "  {:40s} = {} ".format("FEP output file", setup.config["fep"]["output"])
    )
    logger.info(
        "  {:40s} = {} ".format(
            "output frequency", setup.config["fep"]["reportInterval"]
        )
    )

    return lambdaValues, fepAtoms


def vdwAlchemyForAMOEBA(system, fepAtoms):
    """
    Apply alchemical method to AMOEBA VdwForce and modify particle parameters.

    Parameters:
    system (openmm.System): The OpenMM system object.
    fepAtoms (list): A list of atom indices to be modified.

    Returns:
    None
    """
    for f in system.getForces():
        if isinstance(f, mm.AmoebaVdwForce):
            f.setAlchemicalMethod(f.Decouple)
            for i in fepAtoms:
                parm = f.getParticleParameters(i)
                parm[4] = True
                f.setParticleParameters(i, *parm)


def vdwAlchemy(system, fepAtoms):
    """
    Modify the van der Waals (vdW) interactions for specified atoms in the system.

    Parameters:
    - system (openmm.System): The OpenMM system object.
    - fepAtoms (list): A list of atom indices for which the vdW interactions should be modified.

    Returns:
    None
    """
    for f in system.getForces():
        if isinstance(f, mm.CustomNonbondedForce):
            idx = 0
            for idx in range(f.getNumPerParticleParameters()):
                if f.getPerParticleParameterName(idx) == "fepAtom":
                    break

            for iatm in fepAtoms:
                prm = list(f.getParticleParameters(iatm))
                prm[idx] = 0
                f.setParticleParameters(iatm, prm)
