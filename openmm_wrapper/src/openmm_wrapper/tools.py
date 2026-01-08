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
import pprint as pp
import openmm as mm
import openmm.app as app


############################################
def debugDictionary(mydict, name=""):
    logger = logging.getLogger("dynamicEntropy")
    n = 36 - len(name) - 2
    delim = "" + name + ": " + "-" * n + "#"
    logger.debug(delim)
    for k, v in mydict.items():
        if v is not None:
            logger.debug("" + name + ": {}".format(k))
            logger.debug(pp.pformat(v))
    logger.debug(delim)


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
    logger = logging.getLogger("dynamicEntropy")
    logger.debug("-" * 36 + "#")
    logger.debug("Checking available platforms ...")
    numPlatforms = mm.Platform.getNumPlatforms()
    name = ""
    speed0 = 0
    for i in range(numPlatforms):
        platform = mm.Platform.getPlatform(i)
        speed = platform.getSpeed()
        if speed > speed0:
            speed0 = speed
            name = platform.getName()

        logger.debug(
            "  {:33s} = {}".format(platform.getName() + " - relative speed", speed)
        )

    return name


def getNumberOfDegreesOfFreedom(system):
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

    logger = logging.getLogger("dynamicEntropy")
    logger.debug("------------------------------------#")
    logger.debug("NDOF: Number of particles {}".format(numParticles))
    logger.debug("NDOF: Number of barostat dof {}".format(Nbaro))
    logger.debug("NDOF: Number of dof {}".format(ndof))
    logger.debug("------------------------------------#")

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


def get_atomic_coordinates_from_topology(settings, topology, coordinates):
    import openmm.unit as unit
    import openmm_wrapper as my

    types = my.inputToList(settings["species"])
    expr = settings["expr"].replace("X", "x").replace("Y", "y").replace("Z", "z")
    positions_list = []
    for atom in topology.atoms():
        x, y, z = coordinates[atom.index].value_in_unit(unit.nanometer)
        if eval(expr) and atom.name in types:
            atom_data = {
                "atom_name": atom.name,
                "atom_number": atom.index + 1,
                "coord": [x, y, z],
            }
            positions_list.append(atom_data)

    return positions_list


def write_pdb_coordinates(file_path, cell, atoms):
    """
    COLUMNS        DATA  TYPE    FIELD        DEFINITION
    -------------------------------------------------------------------------------------
     1 -  6        Record name   "ATOM  "
     7 - 11        Integer       serial       Atom  serial number.
    13 - 16        Atom          name         Atom name.
    17             Character     altLoc       Alternate location indicator.
    18 - 20        Residue name  resName      Residue name.
    22             Character     chainID      Chain identifier.
    23 - 26        Integer       resSeq       Residue sequence number.
    27             AChar         iCode        Code for insertion of residues.
    31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)     occupancy    Occupancy.
    61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    77 - 78        LString(2)    element      Element symbol, right-justified.
    79 - 80        LString(2)    charge       Charge  on the atom.

    """
    if len(atoms) > 99999:

        def func(x):
            return hex(x + 1)[2:].upper()

    else:

        def func(x):
            return x + 1

    with open(file_path, "w") as pdb_file:
        line = "CRYST1 {} {} {} {} {} {}\n".format(*cell[0:3], *cell[3:6] * 180 / np.pi)
        pdb_file.write(line)
        for i in range(len(atoms)):
            atom_number = func(i)
            atom_type = atoms[i]["atom_name"]
            coord = atoms[i]["coord"]
            line = "{:<6s}{:5} {:<4} {:3}  {:4}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}           {:<2}".format(
                "ATOM",
                atom_number,
                atom_type,
                "UNK",
                atom_number,
                *coord,
                1.00,
                0.00,
                atom_type.upper(),
            )
            pdb_file.write(line)

    return


def getAtoms(config, topology):
    IDlst = []
    if "residue" in config and len(config["residue"]) > 0:
        if isinstance(config["residue"], str):
            config["residue"] = config["residue"].split(",")
        for a in topology.atoms():
            if a.residue.name in config["residue"]:
                IDlst.append(a.index)

    if "species" in config and len(config["species"]) > 0:
        if isinstance(config["species"], str):
            config["species"] = config["species"].split(",")
        for iatm in topology.atoms():
            if iatm.name in config["species"]:
                IDlst.append(iatm.index)

    if "indices" in config:
        if isinstance(config["indices"], list):
            if len(config["indices"]) == 0:
                config["indices"] = []
        elif isinstance(config["indices"], str):
            config["indices"] = [int(x) for x in config["indices"].split(",")]
        elif isinstance(config["indices"], int):
            config["indices"] = [config["indices"]]

        IDlst += config["indices"]

    return IDlst


def createBasicContext(system, setup):
    context = mm.Context(
        system,
        mm.VerletIntegrator(0),
        mm.Platform.getPlatformByName(setup.config["basic"]["Platform"]),
        setup.properties,
    )
    return context


def inputToList(input):
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


def check_required_keywords(my_dict, list_of_keywords, fromfunc=""):
    import sys
    import logging

    logger = logging.getLogger("dynamicEntropy")

    for k in list_of_keywords:
        if k not in my_dict:
            logger.error(f"Missing keyword in function {fromfunc}: {k}")
            sys.error(1)
