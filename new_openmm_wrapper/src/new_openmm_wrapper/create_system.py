### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import logging
import pathlib

import openmm as mm
import openmm.app as app
import new_openmm_wrapper as my


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


def createSystem(setup):
    """
    Create a new system using the forcefield and PDB files

    Parameters
    ----------
        setup:
            object with the configuration parameters
    """
    my.pretty_log()
    # Read the forcefield file(s)
    forcefield = my.readForcefield(setup.config)

    # Open the PDB file and create the openMM modeller object
    my.pretty_log(
        {"Opening coordinates file": setup.config["input"]["coordinates"]}, sep=True
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
        my.pretty_log(
            "New coordinates written to newCoord.pdb", title="Adding virtual site"
        )
        app.PDBFile.writeFile(
            modeller.topology, modeller.positions, open("newCoord.pdb", "w")
        )

    # Delete water
    if setup.config["input"]["dry"] is not None:
        my.pretty_log("Deleting water from the system")
        if setup.config["input"]["dry"] is True:
            modeller.deleteWater()
        else:
            toDelete = [
                r
                for r in modeller.topology.residues()
                if r.name == setup.config["input"]["dry"]
            ]
            modeller.delete(toDelete)

    # Extract the Forcefield settings
    prmFF = extracParametersFF(setup.config["forcefield"])
    my.pretty_log(prmFF, title="Forcefield settings:", sep=True)

    # Check that the residue and atoms names for the water are the ones
    # expected by openMM for rigid water models
    if prmFF["rigidWater"]:
        checkWaterNames(modeller.topology)
    # Ensure the constraints method is correct
    prmFF["constraints"] = fixConstraintsMethod(prmFF["constraints"])

    my.pretty_log("Creating the system object", sep=True)
    system = forcefield.createSystem(modeller.topology, **prmFF)

    if prmFF["rigidWater"] or prmFF["constraints"] is not None:
        my.pretty_log("Checking constraints")
        checkConstraints(system, modeller.topology)

    ##################################################
    my.pretty_log("Creating internal system lists", logger="debug")
    setup.setAtomTypesFromFF(forcefield)
    setup.createListOfAtoms(forcefield, modeller.topology)
    setup.setBondsList(modeller.topology)
    if setup.config["forcefield"]["isAMOEBA"]:
        setup.setExclusionsList(system, nbonds=4)
    else:
        setup.setExclusionsList(system, nbonds=3)

    ##################################################
    # Change mixing rules of oplsAA forcefield is used
    # if setup.config["forcefield"]["isOPLS"]:
    #     my.pretty_log("Creating OPLS mixing rules", logger="debug")
    #     my.OPLS_MixingRules(setup, system)

    # Add restratints if required
    if setup.config["restraint"] is not None:
        my.addRestraints(
            setup.config["restraint"],
            system,
            modeller.topology,
            coordinates=modeller.positions,
        )

    ##################################################
    # global_vars_before_execution = dict(globals())
    if not setup.config["forcefield"]["skipScript"]:
        my.pretty_log("Creating custom forcefields from the XML file")
        for s in forcefield._scripts:
            exec(s, globals())

        for v in ["nonbondedCutoff", "switchDistance"]:
            if v not in globals():
                globals()[v] = setup.config["forcefield"][v]

        if "customForceFieldParameters" in globals():
            # Get customForceFields from globals if it exists, otherwise None
            custom_force_fields = globals().get("customForceFields")

            if custom_force_fields is not None:
                customForceFieldsX = my.createCustomForceFields(
                    customFunctionsFromFile=custom_force_fields
                )
            else:
                customForceFieldsX = my.createCustomForceFields()

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
        my.pretty_log("Custom forcefields have already been created")

    ##################################################

    # if setup.config["forcefield"]["modify"] is not None:
    #     logger.info(
    #         "  {:40s} = {}".format("Modify forcefield", setup.config["forcefield"]["modify"])
    #     )
    #     my.changeAtomsProperties(system, modeller.topology, setup.config["forcefield"]["modify"])

    ##################################################
    if not setup.config["md"]["CMMotionRemover"]:
        cm_remover_type = "disabled"
    else:
        cm_remover_type = setup.config["md"]["CMMotionRemover"].get("type")
        cm_remover_freq = setup.config["md"]["CMMotionRemover"].get("reportInterval")

    # Handle CMMotionRemover - iterate once
    for n, force in enumerate(system.getForces()):
        if force.getName() == "CMMotionRemover":
            if cm_remover_type != "internal":
                system.removeForce(n)
                my.pretty_log("Centre of Mass Motion Remover is Disabled")
            else:
                my.pretty_log(
                    f"Centre of Mass Motion is removed every {force.getFrequency()} steps",
                    logger="warning",
                )
            break

    # # Set force groups
    # log_str = {}
    # for n, f in enumerate(system.getForces()):
    #     f.setForceGroup(n + 1)
    #     log_str[n + 1] = f.getName()
    # my.pretty_log(
    #     title="Creating force groups:", data=log_str, align_width=0, logger="debug"
    # )

    # Compute system properties for basic checks
    my.checkSystem(system, modeller.topology, forcefield, pdb.positions)

    return modeller, system


def checkWaterNames(topology):
    """
    Check that the residue and atoms names for the water are the ones
    expected by openMM for rigid water models

    Parameters
    ----------
        topology:
            openMM topology object
    """
    import sys

    logger = logging.getLogger("dynamicEntropy")
    for res in topology.residues():
        if res.name in ["HOH", "WAT", "TIP3"]:
            expected_atom_names = ("O", "H1", "H2")
            actual_atom_names = tuple(atom.name for atom in res.atoms())
            if actual_atom_names != expected_atom_names:
                logger.error(
                    f"Water residue '{res.name}' has unexpected atom names: {actual_atom_names}.\nExpected names are: {expected_atom_names}."
                )
                sys.exit(1)
    logger.debug("Water residue and atom names are correct for a rigid water model")
    return


def checkConstraints(system, topology):
    """
    Check that the residue and atoms names for the water are the ones
    expected by openMM for rigid water models

    Parameters
    ----------
        topology:
            openMM topology object
    """
    logger = logging.getLogger("dynamicEntropy")

    # Count water molecules
    num_waters = sum(
        1 for res in topology.residues() if res.name in ["HOH", "WAT", "TIP3"]
    )
    expected_constraints = num_waters * 3  # 3 constraints per rigid water

    # Count constraints
    num_constraints = system.getNumConstraints()

    logger.debug(f"Total constraints applied: {num_constraints}")
    logger.debug(f"Number of water molecules: {num_waters}")
    logger.debug(f"Expected water constraints: {expected_constraints}")

    # Check constraint details
    logger.debug("\nFirst 10 constraints:")
    for i in range(min(10, num_constraints)):
        atom1, atom2, distance = system.getConstraintParameters(i)
        logger.debug(f"Constraint {i}: atoms {atom1}-{atom2}, distance {distance}")

    return


def fixConstraintsMethod(str):
    """
    Docstring for fixConstraintsMethod

    :param str: Description
    """
    if str is None:
        return None
    elif str.lower() == "none":
        return None
    elif str.lower() == "allbonds":
        return app.AllBonds
    elif str.lower() == "hbonds":
        return app.HBonds
    elif str.lower() == "hangles":
        return app.HAngles
    elif str in [app.AllBonds, app.HBonds, app.HAngles]:
        return str
    else:
        raise ValueError("Invalid constraints method")
