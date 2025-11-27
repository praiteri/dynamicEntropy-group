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
import pathlib

import openmm as mm
import openmm.app as app
import openmm_wrapper as my


def extracParametersFF(myFF):
    """
    Ensure that the required forcefield parameters are correctly set.
    It removes the AMOEBA specific parameters if the forcefield is not AMOEBA.
    It removes  the ewaldErrorTolerance parameter if the nonbondedMethod is not PME, Ewald or LJPME.

    Parameters
    ----------
        myFF: dictionary with the forcefield parameters
    Returns:
        dictionary with the forcefield parameters that are required to create the system.
    """
    myList = [
        "nonbondedCutoff",
        "ewaldErrorTolerance",
        "nonbondedMethod",
        "useDispersionCorrection",
        "constraints",
        "rigidWater",
    ]
    if myFF["isAMOEBA"]:
        myList.append("polarization")
        myList.append("mutualInducedTargetEpsilon")
    else:
        myList.append("switchDistance")

    if myFF["nonbondedMethod"] in [app.PME, app.Ewald, app.LJPME]:
        assert (
            myFF["ewaldErrorTolerance"] <= 1e-5
        ), f"Ewald error tolerance must be less than 1e-5. Current value is {myFF['ewaldErrorTolerance']}"
    elif myFF["nonbondedMethod"] in [
        app.CutoffPeriodic,
        app.CutoffNonPeriodic,
        app.NoCutoff,
    ]:
        myList.remove("ewaldErrorTolerance")
    else:
        raise Exception("Unknown nonbonded method")

    return {x: myFF[x] for x in myList}


def createSystem(setup, first=True):
    """
    Create a new system using the forcefield and PDB files

    Parameters
    ----------
        setup:
            object with the configuration parameters
        first:
            boolean to print the system info only once
            also used to set the atom types and bonds list...
    """
    logger = logging.getLogger("dynamicEntropy")

    # Read the forcefield file(s)
    forcefield = my.readForcefield(setup.config)

    # Open the PDB file and create the openMM modeller object
    if first:
        logger.critical("Opening coordinates ...")
        logger.info(
            "  {:40s} = {}".format("file", setup.config["input"]["coordinates"])
        )

    coordFile = setup.config["input"]["coordinates"]
    file_extension = pathlib.Path(coordFile).suffix
    if file_extension == ".pdb":
        pdb = app.PDBFile(coordFile)
    elif file_extension == ".pdbx":
        pdb = app.PDBxFile(coordFile)
    else:
        raise Exception("Unsupported coordinates file type - " + coordFile)
    modeller = app.modeller.Modeller(pdb.topology, pdb.positions)

    # Add virtual sites if necessary
    modeller.addExtraParticles(forcefield)
    if len(pdb.positions) != len(modeller.positions):
        app.PDBFile.writeFile(
            modeller.topology, modeller.positions, open("newCoord.pdb", "w")
        )
        if first:
            logger.debug("------------------------------------#")
            logger.debug("NEWSYSTEM: Virtual sites added ...")
            logger.debug("NEWSYSTEM:   New coordinates written to newCoord.pdb")
            logger.debug("------------------------------------#")

    # Delete water
    if setup.config["input"]["dry"] is not None:
        logger.info(
            "  {:40s} = {}".format("Delete water", setup.config["input"]["dry"])
        )
        if setup.config["input"]["dry"] is True:
            modeller.deleteWater()
        else:
            toDelete = [
                r
                for r in modeller.topology.residues()
                if r.name == setup.config["input"]["dry"]
            ]
            modeller.delete(toDelete)

    # Create the system object
    if setup.config["forcefield"]["switchDistance"] is not None:
        sw = float(setup.config["forcefield"]["switchDistance"])
        if sw > 0:
            setup.config["forcefield"]["switchDistance"] = sw
        else:
            setup.config["forcefield"]["switchDistance"] = None

    # Extract the Forcefield settings
    prmFF = extracParametersFF(setup.config["forcefield"])

    # Screen output
    if first:
        setup.dumpParametersFF(prmFF)

    # Create the system
    system = forcefield.createSystem(modeller.topology, **prmFF)

    ##################################################
    if first:
        setup.setAtomTypesFromFF(forcefield)
        setup.createListOfAtoms(forcefield, modeller.topology)
        setup.setBondsList(modeller.topology)
        if setup.config["forcefield"]["isAMOEBA"]:
            setup.setExclusionsList(system, nbonds=4)
        else:
            setup.setExclusionsList(system, nbonds=3)
        # setup.setExclusionsList(system, nbonds=3)

    ##################################################
    # Change mixing rules of oplsAA forcefield is used
    if setup.config["forcefield"]["isOPLS"]:
        my.OPLS_MixingRules(setup, system)

    # Add restratints if required
    if setup.config["restraint"] is not None:
        my.addRestraints(setup.config["restraint"], system, modeller.topology)

    ##################################################
    # global_vars_before_execution = dict(globals())
    if not setup.config["forcefield"]["skipScript"]:
        logger.debug("Creating custom forcefields from the XML file ...")
        for s in forcefield._scripts:
            exec(s, globals())

        for v in ["nonbondedCutoff", "switchDistance"]:
            if not v in globals():
                globals()[v] = setup.config["forcefield"][v]

        if "customForceFieldParameters" in globals():
            if "customForceFields" not in globals():
                customForceFieldsX = my.customForceFields()
            else:
                customForceFieldsX = my.customForceFields(
                    customFunctionsFromFile=customForceFields
                )

            my.customNonbondedForce(
                system,
                setup,
                customForceFieldParameters,
                customForceFieldsX,
                nonbondedCutoff,
            )

            my.customNonbonded14Force2(
                system,
                setup,
                customForceFieldParameters,
                customForceFieldsX,
                setup.config["forcefield"]["lj14Scale"],
            )

        if "customBondForceParameters" in globals():
            my.customBondForce(system, setup, customBondForceParameters)

        if "customImproperParameters" in globals():
            my.customDistanceImproperForce(
                system, setup, forcefield, modeller.topology, customImproperParameters
            )

        if "customCosineAngle" in globals():
            my.customCosineAngleForce(
                system, setup, forcefield, modeller.topology, customCosineAngle
            )
    else:
        logger.debug("Custom forcefields have already been created ...")

    ##################################################

    if setup.config["forcefield"]["modify"] is not None:
        logger.info(
            "  {:40s} = {}".format(
                "Modify forcefield", setup.config["forcefield"]["modify"]
            )
        )
        my.changeAtomsProperties(
            system, modeller.topology, setup.config["forcefield"]["modify"]
        )

    ##################################################

    if setup.config["forcefield"]["CMMotionRemover"] is False:
        for n, force in enumerate(system.getForces()):
            if force.getName() == "CMMotionRemover":
                system.removeForce(n)
                break
        logger.warning("Centre of Mass Motion Remover is Disabled")
    else:
        logger.warning("Centre of Mass Motion Remover is Enabled")

    if first:
        logger.debug("------------------------------------#")
        logger.debug("NEWSYSTEM: Adding force groups ...")
        for n, f in enumerate(system.getForces()):
            logger.debug("NEWSYSTEM:   {} -> {}".format(n + 1, f.getName()))
            f.setForceGroup(n + 1)
        logger.debug("------------------------------------#")
        # Compute system properties for basic checks
        my.checkSystem(system, modeller.topology, forcefield, pdb.positions)
    else:
        for n, f in enumerate(system.getForces()):
            f.setForceGroup(n + 1)

    # Remove dispersion correction
    # if not setup.config["forcefield"]["useDispersionCorrection"]:
    #     for force in system.getForces():
    #         if force.getName() == "NonbondedForce":
    #             force.setUseDispersionCorrection(False)

    return modeller, system
