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
import openmm.unit as unit
import numpy as np


def customNonbonded14Force2(
    system, setup, interactionsList, customForceFields, lj14scale
):
    logger = logging.getLogger("dynamicEntropy")

    if lj14scale <= 0.0:
        return
    logger.debug("FF_14: lj14scale = {}".format(lj14scale))

    # Get nonBondedMethod for the interactions
    # 0 -> NoCutoff; 4 -> PME
    for n, f in enumerate(system.getForces()):
        if isinstance(f, mm.NonbondedForce) or isinstance(f, mm.AmoebaMultipoleForce):
            nbMethod = f.getNonbondedMethod()

    # list of atom types
    listOfAtoms = setup.getListOfAtoms()

    # List of bonds in the system
    bonds = setup.getBondsList()

    # List of 1-4 interactions in the system
    exclusions = app.forcefield._findExclusions(bonds, 3, system.getNumParticles())
    list14 = [x for x in exclusions if x[2] == 3]
    if len(list14) == 0:
        return

    # Loop over all th custom forcefields defined
    for cft in customForceFields:
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
                if cft["type"] != iList["ff"]:
                    continue
                if (iList["type1"] == a1 and iList["type2"] == a2) or (
                    iList["type1"] == a2 and iList["type2"] == a1
                ):
                    p = [iList["params"][x] for x in cft["params"]]
                    force14.addBond(i1, i2, p)
                    numInt += 1

        if numInt > 0:
            force14.setName(cft["type"] + "_14_vdwl")
            system.addForce(force14)
        else:
            del force14


# def customNonbonded14Force(sys,topology,setup,interactionsList,lj14scale):
#     """
#     Add 1-4 interactions for custom LJ non bonded forces that are used
#     to override the default mixing rules.
#     Parameters
#     ----------
#         lj14scale: scaling factor for the LJ interactions
#     """
#     # list of atom types
#     listOfAtoms = setup.getListOfAtoms()

#     # List of bonds in the system
#     bonds = setup.getBondsList()

#     # List of 1-4 interactions in the system
#     exclusions = app.forcefield._findExclusions(bonds,3,sys.getNumParticles())
#     list14 = [x for x in exclusions if x[2] == 3]

#     # Custom 1-4 interactions
#     expression = "{}*e*((s/r)^12-(s/r)^6); ".format(4*lj14scale)
#     lj14 = mm.CustomBondForce(expression)
#     lj14.addPerBondParameter('s')
#     lj14.addPerBondParameter('e')

#     numInt = 0
#     for i in list14:
#         i1 = i[0]
#         i2 = i[1]
#         a1 = listOfAtoms[i1]["type"]
#         a2 = listOfAtoms[i2]["type"]

#         for iList in interactionsList:
#             if "lj" not in iList["ff"]: continue
#             if (iList["type1"]==a1 and iList["type2"]==a2) or (iList["type1"]==a2 and iList["type2"]==a1):
#                 eps = iList["params"]['e']
#                 sig = iList["params"]['s']
#                 if "conv" in iList:
#                     raise Exception("conv has been removed, add conversion to the parameters")
#                     # eps = eps*iList["conv"]['e']
#                     # sig = sig*iList["conv"]['s']
#                 lj14.addBond(i1, i2, [sig, eps])
#                 numInt += 1

#     if numInt > 0:
#         lj14.setName("lj14")
#         sys.addForce(lj14)
#     else:
#         del lj14
