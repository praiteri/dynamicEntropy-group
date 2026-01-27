### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

def input_to_list(input):
    if type(input) is list:
        return input

    input = str(input)
    result = ""
    for char in input:
        if char not in "[](){}":
            result += char
    return result.split(",")


def replace_entries_recursive(dict1, dict2, main=True):
    import openmm.unit as unit
    import openmm.app as app
    import ast

    for key1, value1 in dict1.items():
        if value1 is None and main:
            continue
        if key1 in dict2:
            if isinstance(value1, dict) and isinstance(dict2[key1], dict):
                replace_entries_recursive(value1, dict2[key1], main=False)
            else:
                if unit.is_quantity(dict2[key1]):
                    dict2[key1] = float(dict1[key1]) * dict2[key1].unit

                elif isinstance(dict2[key1], type(app.forcefield.PME)):
                    dict2[key1] = eval("app." + dict1[key1])

                else:
                    try:
                        dict2[key1] = ast.literal_eval(dict1[key1])
                    except:
                        dict2[key1] = dict1[key1]
        else:
            raise ValueError(f"Key '{key1}' not found in {dict2}")


def locate_and_modify_key(data, target_key, new_value):
    import ast

    target_key = target_key.strip()
    if isinstance(data, dict):
        for key, value in data.items():
            if key == target_key:
                data[key] = modifyParameter(data[key], new_value)
                return True
            elif isinstance(value, dict):
                locate_and_modify_key(value, target_key, new_value)
        data[target_key] = ast.literal_eval(new_value)
    return False


def modifyParameter(key, value):
    import ast
    import openmm as mm
    import openmm.app as app

    if type(key) is str or key is None:
        if value == "None":
            key = None
        else:
            key = value
    elif type(key) is mm.unit.quantity.Quantity:
        key = float(value) * key.unit
    elif isinstance(key, type(app.forcefield.PME)):
        key = eval("app." + value)
    elif type(key) is list:
        value = "".join(char for char in value if char not in "[]")
        key = value.split(",")
    else:
        key = ast.literal_eval(value)
    return key


def checkPlatforms():
    """
    Check the available Platforms are return the name of the faster one.
    """
    import openmm as mm
    import new_openmm_wrapper as my

    numPlatforms = mm.Platform.getNumPlatforms()
    name = ""
    speed0 = 0

    dbg = {}
    for i in range(numPlatforms):
        platform = mm.Platform.getPlatform(i)
        speed = platform.getSpeed()
        if speed > speed0:
            speed0 = speed
            name = platform.getName()

        dbg[platform.getName()] = speed
    my.pretty_log(
        dbg,
        logger="DEBUG",
        title="Checking available platforms and relative speeds:",
        sep=True,
    )

    return name


def get_number_of_degrees_of_freedom(system, logger=None):
    import new_openmm_wrapper as my
    import openmm as mm

    # Number of particles in the system
    numParticles = system.getNumParticles()

    # Remove any virtual sites from the particles' count
    for i in range(0, numParticles):
        if system.isVirtualSite(i):
            numParticles -= 1

    # Number of constraints in the system
    Nconst = system.getNumConstraints()

    # Number of degrees of freedom for the barostat
    Nbaro = 0
    for force in system.getForces():
        if isinstance(force, mm.MonteCarloBarostat):
            Nbaro = 1
        elif isinstance(force, mm.MonteCarloAnisotropicBarostat):
            if force.getScaleX():
                Nbaro += 1
            if force.getScaleY():
                Nbaro += 1
            if force.getScaleZ():
                Nbaro += 1
        elif isinstance(force, mm.MonteCarloFlexibleBarostat):
            Nbaro = 6

    ndof = numParticles * 3 - 3 - Nconst + Nbaro

    if logger is None:
        return ndof

    log_str = {}
    log_str["Number of particles"] = numParticles
    log_str["Number of constraints"] = Nconst
    log_str["Number of barostat dof"] = Nbaro
    log_str["Total degrees of freedom"] = ndof
    my.pretty_log(log_str, title="Degree of freedom analysis:", logger=logger)

    return ndof


def read_pdb_coordinates(file_path, unit="A"):
    import numpy as np

    """
    Read atomic coordinates from a PDB file and return them as a list of dictionaries.
    Each dictionary contains the following keys:
    - 'atom_name': Atom name (e.g., 'CA', 'C', 'N')
    - 'coord': X,Y,Z coordinates
    """
    coordinates = []

    if unit.lower() == "nm":

        def func(x):
            return x / 10

    else:

        def func(x):
            return x

    cell = None
    with open(file_path, "r") as pdb_file:
        for line in pdb_file:
            if line.startswith("CRYST"):
                val = line.split()[1:7]
                cell = np.array([float(x) for x in val])
                cell[0:3] = func(cell[0:3])
                cell[3:6] = cell[3:6] * np.pi / 180.0

            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    atom_number = int(line[6:11].strip())
                except:
                    atom_number += 1
                atom_name = line[12:16].strip()
                # record_name = line[0:6].strip()
                # residue_name = line[17:20].strip()
                # chain_id = line[21:22].strip()
                # residue_number = int(line[22:26].strip())
                x = func(float(line[30:38].strip()))
                y = func(float(line[38:46].strip()))
                z = func(float(line[46:54].strip()))

                atom_data = {
                    "atom_name": atom_name,
                    "atom_number": atom_number,
                    "coord": [x, y, z],
                }

                coordinates.append(atom_data)

    return coordinates, cell


def safe_eval_expr(expr, x, y, z):
    """
    Safely evaluate a spatial expression with x, y, z coordinates.

    Args:
        expr: String expression like "x > 2.0 and y < 5.0"
        x, y, z: Coordinate values

    Returns:
        Boolean result of the expression
    """
    from simpleeval import simple_eval

    return simple_eval(expr, names={"x": x, "y": y, "z": z})


def get_atomic_coordinates(atom_list, coordinates):
    import openmm.unit as unit
    import new_openmm_wrapper as my

    positions_list = []
    for atom in atom_list:
        x, y, z = coordinates[atom.index].value_in_unit(unit.nanometer)
        atom_data = {
            "atom_name": atom.name,
            "atom_number": atom.index + 1,
            "coord": [x, y, z],
        }
        positions_list.append(atom_data)

    return positions_list


# def getAtoms(config, topology):
#     IDlst = []
#     if "residue" in config and len(config["residue"]) > 0:
#         if isinstance(config["residue"], str):
#             config["residue"] = config["residue"].split(",")
#         for a in topology.atoms():
#             if a.residue.name in config["residue"]:
#                 IDlst.append(a.index)

#     if "species" in config and len(config["species"]) > 0:
#         if isinstance(config["species"], str):
#             config["species"] = config["species"].split(",")
#         for iatm in topology.atoms():
#             if iatm.name in config["species"]:
#                 IDlst.append(iatm.index)

#     if "indices" in config:
#         if isinstance(config["indices"], list):
#             if len(config["indices"]) == 0:
#                 config["indices"] = []
#         elif isinstance(config["indices"], str):
#             config["indices"] = [int(x) for x in config["indices"].split(",")]
#         elif isinstance(config["indices"], int):
#             config["indices"] = [config["indices"]]

#         IDlst += config["indices"]

#     return IDlst


def get_atoms_from_topology(config, topology, coordinates=None):
    import sys
    import new_openmm_wrapper as my
    import openmm.unit as unit

    selection_dict = config.get("select")

    if selection_dict is None:
        my.pretty_log(config, title="Cannot select atoms from command", logger="error")
        sys.exit(1)

    IDlst = []
    for ilk, value in selection_dict.items():
        if ilk in ["residue", "resname"]:
            if isinstance(value, str):
                selections = value.split(",")
            elif isinstance(value, list):
                selections = value
            else:
                my.pretty_log(f"Cannot understand residue selection ({value})", logger="error")
                sys.exit(1)
            lsel = [iatm.index for iatm in topology.atoms() if iatm.residue.name in selections]
            IDlst += lsel

        elif ilk in "species":
            if isinstance(value, str):
                selections = value.split(",")
            elif isinstance(value, list):
                selections = value
            else:
                my.pretty_log(f"Cannot understand species selection {value}", logger="error")
                sys.exit(1)
            lsel = [iatm.index for iatm in topology.atoms() if iatm.name in selections]
            IDlst += lsel

        elif ilk in ["index", "indices"]:
            if isinstance(value, str):
                selections = [int(x) for x in value.split(",")]
            elif isinstance(value, list):
                selections = [int(x) for x in value]
            elif isinstance(value, int):
                selections = [value]
            else:
                my.pretty_log(f"Cannot understand index selection {value}", logger="error")
                sys.exit(1)
            IDlst += selections

        elif ilk in "positions":
            if coordinates is None:
                my.pretty_log("Missing coordinates, cannot do position selection", logger="error")
                sys.exit(1)
        else:
            my.pretty_log(f"Cannot understand selection {ilk}", logger="error")
            sys.exit(1)

    for ilk, value in selection_dict.items():
        if ilk in "positions":
            expr = value.replace("X", "x").replace("Y", "y").replace("Z", "z")
            if len(IDlst) == 0:
                IDlst = [atom.index for atom in topology.atoms()]

            new_IDlst = []
            for iatm in IDlst:
                x, y, z = coordinates[iatm].value_in_unit(unit.nanometer)
                if safe_eval_expr(expr, x, y, z):
                    new_IDlst.append(iatm)

    return IDlst


def check_required_keywords(my_dict, list_of_keywords, fromfunc=""):
    import sys
    import logging

    logger = logging.getLogger("dynamicEntropy")

    for k in list_of_keywords:
        if k not in my_dict:
            logger.error(f"Missing keyword in function {fromfunc}: {k}")
            sys.exit(1)


def is_debug(config):
    if "debug" in config and config["debug"]:
        return True
    return False
