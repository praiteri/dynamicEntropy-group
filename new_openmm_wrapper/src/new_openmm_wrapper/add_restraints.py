### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import ast

import openmm as mm
import new_openmm_wrapper as my


def make_list_IDs(IDlst, topology):
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


def positionsRestraint(forceName, settings, system, topology, coordinates=None):
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
    import sys

    my.pretty_log(f"{forceName}:", indent=1)

    variables = my.input_to_list(settings["par"])
    values = [float(x) for x in my.input_to_list(settings["val"])]

    globalParameters = {}
    if "global" in settings:
        for p in settings["global"].split(","):
            idx = variables.index(p)
            globalParameters[variables[idx]] = values[idx]
            variables.pop(idx)
            values.pop(idx)

    # if "fullexp" in settings:
    #     expression = settings["fullexp"]

    # elif "exp" in settings:
    #     energy = settings["exp"]
    #     distance = " d = periodicdistance(x,y,z,x0,y0,z0)"

    # else:
    #     energy = "k*d^2"
    #     distance = " d = periodicdistance(x,y,z,x0,y0,z0)"
    #     f = "global" in settings and "d0" in settings["global"]
    #     if "d0" in variables or f:
    #         energy = energy.replace("d", "(max(0,d-d0))")

    # This is not a per-atom positional restraint
    if "file" not in settings and "initial" not in settings:
        atomsList = my.get_atoms_from_topology(settings, topology, coordinates)
        parameters = [values] * len(atomsList)

    # This is a positional restraint
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

        if "file" in settings and settings["file"] is not None:
            try:
                atomic_coordinates, _ = my.read_pdb_coordinates(
                    settings["file"], unit="nm"
                )
            except FileNotFoundError:
                raise Exception("Restraint position files missing - ", settings["file"])
            except Exception as e:
                raise Exception("Error reading restraint position file - ", e)

        elif coordinates is not None and "select" in settings:
            atomsList = my.get_atoms_from_topology(settings, topology, coordinates)

            atomic_coordinates = my.get_atomic_coordinates(atomsList, coordinates)
            my.pretty_log(
                f"  {len(atomic_coordinates)} selected atoms from topology ..."
            )
        else:
            sys.exit(1)

        for atom in atomic_coordinates:
            atomsList.append(atom["atom_number"] - 1)
            parameters.append(values + atom["coord"])
        assert len(atomsList) > 0, "No atoms in the restraint file ({})".format(
            settings["file"]
        )

    ##############

    if "fullexp" in settings:
        expression = settings["fullexp"]
    else:
        # Default energy expression
        energy = settings.get("exp", "k*d^2")

        # Apply d0 offset if available
        if "d0" in variables or ("global" in settings and "d0" in settings["global"]):
            energy = energy.replace("d", "(max(0,d-d0))")

        # Build distance calculation
        distance = "d = periodicdistance(x,y,z,x0,y0,z0)"

        # Replace unused reference points with their current values
        if "x0" not in variables:
            distance = distance.replace("x0", "x")
        if "y0" not in variables:
            distance = distance.replace("y0", "y")
        if "z0" not in variables:
            distance = distance.replace("z0", "z")

        expression = f"{energy} ; {distance}"

    # if "ti" in settings and settings["ti"]:
    #     if "lambdaTI" not in expression:
    #         expression = "lambdaTI*" + expression

    # if "fep" in settings and settings["fep"]:
    #     if "lambdaFEP" not in expression:
    #         expression = "lambdaFEP*" + expression

    force = mm.CustomExternalForce(expression)

    if len(globalParameters) > 0:
        for k, v in globalParameters.items():
            force.addGlobalParameter(k, v)

    # if "lambdaTI" in expression:
    #     force.addGlobalParameter("lambdaTI", 0.0)

    if "lambdaFEP" in expression:
        force.addGlobalParameter("lambdaFEP", 0.0)
    if "lambda" in expression:
        force.addGlobalParameter("lambda", 0.0)

    for x in variables:
        force.addPerParticleParameter(x)

    for i, j in zip(atomsList, parameters):
        if isinstance(j, tuple):
            j = list(j)
        elif not isinstance(j, list):
            j = [j]
        force.addParticle(i, j)

    force.setName(forceName)

    # n = system.getNumForces()
    # force.setForceGroup(n + 1)
    system.addForce(force)

    ######
    my.pretty_log({"Number of atoms to be restrained": len(atomsList)}, indent=2)
    my.pretty_log("Expression", indent=2)
    for x in expression.split(";"):
        my.pretty_log(x.strip(), indent=3)
    my.pretty_log({"Input variables": settings["par"]}, indent=2)
    my.pretty_log({"Input values": settings["val"]}, indent=2)
    my.pretty_log({"Local variables": variables}, indent=2)
    if len(globalParameters) > 0:
        my.pretty_log({"Global variables": list(globalParameters.keys())}, indent=2)

    if "lambdaTI" in expression:
        my.pretty_log("Restraint can be used for TI calculations", indent=2)

    if "lambdaFEP" in expression:
        my.pretty_log("Restraint can be used for FEP calculations", indent=2)


def addRestraints(cmd, system, topology, coordinates=None):
    """
    Adds restraints to the system based on the given command and configuration.

    Parameters:
        - cmd (dict or list): The command or configuration for adding restraints.
        - system: The OpenMM system object.
        - topology: The OpenMM topology object.

    Returns:
        None
    """

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

    my.pretty_log(title="Adding restraints", sep=True)

    for name, config in cmd.items():
        positionsRestraint(name, config, system, topology, coordinates)
