### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import logging
import pprint as pp
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np


def customNonbonded14Force2(system, setup, interactionsList, customForceFields, lj14scale):
    logger = logging.getLogger("dynamicEntropy")

    if lj14scale <= 0.0:
        return
    logger.debug("FF_14: lj14scale = {}".format(lj14scale))

    # Get nonBondedMethod for the interactions
    # 0 -> NoCutoff; 4 -> PME
    # for n, f in enumerate(system.getForces()):
    #     if isinstance(f, mm.NonbondedForce) or isinstance(f, mm.AmoebaMultipoleForce):
    #         nbMethod = f.getNonbondedMethod()

    # list of atom types
    listOfAtoms = setup.getListOfAtoms()

    # List of bonds in the system
    bonds = setup.getBondsList()

    # List of 1-4 interactions in the system
    exclusions = app.forcefield._findExclusions(bonds, 3, system.getNumParticles())
    list14 = [x for x in exclusions if x[2] == 3]
    if len(list14) == 0:
        return

    # Loop over all the custom forcefields defined
    for ff_type, cft in customForceFields.items():
        func = cft["expression"].split(";")
        expression = "{}*({})".format(lj14scale, func[0])

        logger.debug("FF_14: expression:")
        logger.debug("FF_14:   {}".format(expression))

        force14 = mm.CustomBondForce(expression)
        for p in cft["params"]:
            force14.addPerBondParameter(p)

        numInt = 0
        for i in list14:
            i1, i2 = i[0:2]
            a1 = listOfAtoms[i1]["type"]
            a2 = listOfAtoms[i2]["type"]

            for iList in interactionsList:
                if ff_type != iList["ff"]:
                    continue
                if (iList["type1"] == a1 and iList["type2"] == a2) or (
                    iList["type1"] == a2 and iList["type2"] == a1
                ):
                    p = [iList["params"][x] for x in cft["params"]]
                    force14.addBond(i1, i2, p)
                    numInt += 1

        if numInt > 0:
            force14.setName(ff_type + "_14_vdwl")
            system.addForce(force14)
        else:
            del force14
