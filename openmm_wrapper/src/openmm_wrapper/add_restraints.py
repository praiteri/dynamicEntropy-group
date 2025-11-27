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
from pprint import *

import openmm as mm
import openmm_wrapper as my
import ast


def makeListIDs(IDlst, topology):
    """
    Create a list of atom indices based on the given atom names.

    Parameters:
        IDlst (list): A list of atom names.
        topology (Topology): The OpenMM Topology object.

    Returns:
        list: A list of atom indices corresponding to the given atom names.
    """
    selection = [atom.index for atom in topology.atoms() if atom.name in IDlst]
    return selection


def positionsRestraint(name, settings, system, topology):
    """
    Apply position restraints to a system based on the given settings.

    Args:
        name (str): The name of the restraint force.
        settings (dict): A dictionary containing the settings for the restraints.
        system (mm.System): The OpenMM System object to which the restraints will be added.
        topology (mm.app.Topology): The OpenMM Topology object representing the system.

    Raises:
        Exception: If the restraint position files are missing.

    Returns:
        None
    """
    logger = logging.getLogger("dynamicEntropy")

    forceName = name
    logger.info("  restraint name: {}".format(forceName))

    variables = my.inputToList(settings["par"])
    values = [float(x) for x in my.inputToList(settings["val"])]

    globalParameters = []
    if "global" in settings:
        for p in settings["global"].split(","):
            idx = variables.index(p)
            globalParameters.append([variables[idx], values[idx]])
            variables.pop(idx)
            values.pop(idx)

    if "fullexp" in settings:
        expression = settings["fullexp"]

    elif "exp" in settings:
        energy = settings["exp"]
        distance = " d = periodicdistance(x,y,z,x0,y0,z0)"

    else:
        energy = "k*d^2"
        distance = " d = periodicdistance(x,y,z,x0,y0,z0)"
        f = "global" in settings and "d0" in settings["global"]
        if "d0" in variables or f:
            energy = energy.replace("d", "(max(0,d-d0))")

    if "file" not in settings:
        atomsList = my.getAtoms(settings, topology)
        parameters = [values] * len(atomsList)
    else:
        atomsList = []
        parameters = []
        # Move the x, y and z at the end of the list
        if "x0" in variables:
            variables.remove("x0")
        if "y0" in variables:
            variables.remove("y0")
        if "z0" in variables:
            variables.remove("z0")
        variables += ["x0", "y0", "z0"]

        try:
            atomic_coordinates, _ = my.read_pdb_coordinates(settings["file"], unit="nm")
        except:
            raise Exception("Restraint position files missing - ", settings["file"])

        for atom in atomic_coordinates:
            atomsList.append(atom["atom_number"] - 1)
            parameters.append(values + atom["coord"])
        assert len(atomsList) > 0, "No atoms in the restraint file ({})".format(
            settings["file"]
        )

    if "fullexp" not in settings:
        if "x0" not in variables:
            distance = distance.replace("x0", "x")
        if "y0" not in variables:
            distance = distance.replace("y0", "y")
        if "z0" not in variables:
            distance = distance.replace("z0", "z")
        expression = energy + " ; " + distance

    if "ti" in settings and settings["ti"]:
        if "lambdaTI" not in expression:
            expression = "lambdaTI*" + expression

    if "fep" in settings and settings["fep"]:
        if "lambdaFEP" not in expression:
            expression = "lambdaFEP*" + expression

    logger.info("  expression: {}".format(expression.replace(";", "\n  ")))
    for x, y in zip(variables, values):
        logger.info("    {} = {}".format(x, y))

    if len(globalParameters) > 0:
        for x in globalParameters:
            logger.info("    {} = {}".format(*x))

    if "file" in settings:
        logger.info("  atoms selected from file: {}".format(settings["file"]))
    elif "indices" in settings:
        logger.info("  selected atoms: {}".format(settings["indices"]))
    elif "species" in settings:
        logger.info("  selected species: {}".format(settings["species"]))

    force = mm.CustomExternalForce(expression)

    if len(globalParameters) > 0:
        for x in globalParameters:
            force.addGlobalParameter(*x)

    if "lambdaTI" in expression:
        force.addGlobalParameter("lambdaTI", 0.0)
        logger.info("  {} can be used for TI calculations".format(forceName))

    if "lambdaFEP" in expression:
        force.addGlobalParameter("lambdaFEP", 0.0)
        logger.info("  {} can be used for FEP calculations".format(forceName))

    for x in variables:
        force.addPerParticleParameter(x)

    logger.info("  Number of atoms selected: {}".format(len(atomsList)))

    for i, j in zip(atomsList, parameters):
        if type(j) is tuple:
            j = [
                *j,
            ]
        elif type(j) is not list:
            j = [j]
        force.addParticle(i, j)

    force.setName(forceName)

    n = system.getNumForces()
    force.setForceGroup(n + 1)
    system.addForce(force)


def addRestraints(cmd, system, topology):
    """
    Adds restraints to the system based on the given command and configuration.

    Parameters:
        - cmd (dict or list): The command or configuration for adding restraints.
        - system: The OpenMM system object.
        - topology: The OpenMM topology object.

    Returns:
        None
    """
    logger = logging.getLogger("dynamicEntropy")
    logger.critical("Adding restraints ...")

    # This to cope with the restraints read from the command line
    # instead of the yaml file
    if isinstance(cmd, list):
        flat_list = [item for sublist in cmd for item in sublist]
        string = ",".join(flat_list)

        index = string.rfind("},{")
        while index != -1:
            string = string[:index] + "," + string[index + 3 :]
            index = string.rfind("},{")
        cmd = ast.literal_eval(string)

    logger.debug("RESTRAINTS: ------------------------#")
    for name, config in cmd.items():
        logger.debug("RESTRAINTS: {}".format(name))
        logger.debug(pformat(config))

        positionsRestraint(name, config, system, topology)
    logger.debug("RESTRAINTS: ------------------------#")


# Your code goes here
