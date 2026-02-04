### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import logging
import openmm as mm
import numpy as np
import re


def customNonbondedForce(
    system, setup, interactionsList, customForceFields, cutoffDistance
):
    """
    Add custom non bonded forces defined in the focefield XML file.

    Parameters
    ----------
    system: System
        the System to simulate

    setup: dict
        simulation setup parameters

    interactionsList:
        list of forcefield pair interactterms defined in the xml file.
        Each entry contains a dictionary for the interactions.
        The parameters and conversion factors must be in the same order
        as the parameters defined in the customForceField entries.
            myFF.append({ "ff":"lj"   , "type1":"OW" , "type2":"OW" , "params":[0.00674  , 3.165492 ] , "conv":[ev2kj , a2nm] })
            myFF.append({ "ff":"lj"   , "type1":"Mg" , "type2":"OW" , "params":[0.001137 , 2.720000 ] , "conv":[ev2kj , a2nm] })

    customForceFields:
        dict of custom defined potentials in the forcefield xml file.

    cutoffDistance:
        cutoff distance for the interactions
    """
    logger = logging.getLogger("dynamicEntropy")

    # Get nonBondedMethod for the interactions
    # 0 -> NoCutoff
    # 4 -> PME

    for n, f in enumerate(system.getForces()):
        if isinstance(f, mm.NonbondedForce) or isinstance(f, mm.AmoebaMultipoleForce):
            nbMethod = f.getNonbondedMethod()

    listOfAtomTypesPresent = setup.getListOfAtomTypesPresent()
    numTypes = len(listOfAtomTypesPresent)
    listOfAtoms = setup.getListOfAtoms()

    logger.debug("Atom types in pdb file")
    logger.debug(" {}".format(listOfAtomTypesPresent))

    # Loop over all the custom forcefields defined
    for ff_type, cft in customForceFields.items():
        # Check if this forcefield type is used
        isRequired = False
        for iList in interactionsList:
            if iList["ff"] == ff_type:
                # isRequired = True
                if (
                    iList["type1"] in listOfAtomTypesPresent
                    and iList["type2"] in listOfAtomTypesPresent
                ):
                    isRequired = True
        if not isRequired:
            continue

        # Skip if the entry has no "expression"
        if "expression" not in cft:
            continue

        logger.debug("FF: Custom forcefield type = {}".format(ff_type))
        # Create a new force object
        expression = cft["expression"]

        logger.debug("FF: expression:")
        for x in expression.split(";"):
            logger.debug("FF:   {}".format(x.replace(" ", "")))
        newForce = mm.CustomNonbondedForce(expression)

        if nbMethod == 0:
            newForce.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        else:
            newForce.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
            newForce.setCutoffDistance(cutoffDistance)

        switchDistance = -1
        if setup.config["forcefield"]["customSwitchDistance"] is not None:
            switchDistance = float(setup.config["forcefield"]["customSwitchDistance"])

        if switchDistance <= 0:
            newForce.setUseSwitchingFunction(False)
        else:
            newForce.setUseSwitchingFunction(True)
        newForce.setSwitchingDistance(switchDistance)

        # Only one paramater per particle is required since the actual parameters are given as matrices
        _ = newForce.addPerParticleParameter("type")

        for a in listOfAtoms:
            newForce.addParticle([a["itype"]])

        # Add the parameters
        for prm in cft["params"]:
            if cft["init"][prm] == 1:
                paramMatrix = np.ones([numTypes, numTypes])
            else:
                paramMatrix = np.zeros([numTypes, numTypes])

            for iList in interactionsList:
                if iList["ff"] == ff_type:
                    if (iList["type1"] in listOfAtomTypesPresent) and (
                        iList["type2"] in listOfAtomTypesPresent
                    ):
                        i = listOfAtomTypesPresent.index(iList["type1"])
                        j = listOfAtomTypesPresent.index(iList["type2"])

                        if prm in iList["params"]:
                            p = iList["params"][prm]
                        else:
                            if "fepprm" in iList:
                                # if prm in iList["fepprm"]:
                                p = iList["fepprm"][prm]
                            else:
                                s = prm[:-1]
                                p = iList["params"][s]

                        if "conv" in iList:
                            raise Exception(
                                "conv has been removed, add conversion to the parameters"
                            )
                            # p *= iList["conv"][prm]
                        paramMatrix[i, j] = paramMatrix[j, i] = p

            debugCustomForcefieldInteractions(prm, listOfAtomTypesPresent, paramMatrix)

            listOfParameters = (paramMatrix).ravel().tolist()
            newForce.addTabulatedFunction(
                prm, mm.Discrete2DFunction(numTypes, numTypes, listOfParameters)
            )

        # Set the exclustions as in the nonBondedForce
        for l in setup.getExclusionsList():
            newForce.addExclusion(*l[0:2])

        # Add the force to the systemtem
        newForce.setName(ff_type + "_vdwl")
        system.addForce(newForce)


def debugCustomForcefieldInteractions(fparam, listOfAtomTypes, paramMatrix):
    import logging

    logger = logging.getLogger("dynamicEntropy")
    logger.debug("FF: parameter = {}".format(fparam))
    str_len = 8
    line = " " * (str_len + 1)
    for xx in range(len(paramMatrix)):
        line += f"{listOfAtomTypes[xx]:{str_len}s} "
    logger.debug(line)
    for xx in range(len(paramMatrix)):
        line = f"{listOfAtomTypes[xx]:{str_len}s} "
        for yy in range(len(paramMatrix)):
            x = paramMatrix[xx, yy]
            if x > 0:
                n = int(np.log10(x))
                m = min(6, 6 - n)
                if n > 6:
                    string = f"{round(x, m):<.{str_len}g}"
                else:
                    string = "{:<.{width}f}".format(round(x, m), width=m)
            else:
                m = 6
                string = "{:.6f}".format(round(x, m))
            if len(string) == 7:
                string += "."
            line += string + " "
        logger.debug(line)
    line = "_" * (len(paramMatrix) * 8 + 7)
    logger.debug(line)
    return
