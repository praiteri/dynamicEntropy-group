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
import numpy as np
import pprint as pp
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
        list of custom defined potentials in the forcefield xml file.
        Init defined the initialisation of the parameters, useful to avoid
        computing 1/0 while doing FEP calculations.
            customForceFields.append({"type":"lj",
                                    "params":["e","s"],
                                    "init":{'e':0,'s':1},
                                    "expression":"4*e*( (s/r)^12 - (s/r)^6 ); \
                                                    e=e(type1, type2); \
                                                    s=s(type1, type2)",
                                    "fep" : "lambdaVdw^2*4*e*(1./(alpha+(r/s)^6)^2-1/(alpha+(r/s)^6))"
                                    })

            customForceFields.append({"type":"AB",
                                    "params":["A","B"],
                                    "init":{'A':0,'B':0},
                                    "expression":"A/d^12 - B/d^6; \
                                                    A=A(type1, type2); \
                                                    B=B(type1, type2)",
                                    "fep" : "lambdaVdw^2*(A*exp(-r/rho) - C/(alpha+r^2)^3)"
                                    })

            customForceFields.append({"type":"buck",
                                    "params":["A","rho","C"],
                                    "init":{'A':0,'rho':1,'C':0},
                                    "expression":"A*exp(-r/rho) - C/r^6; \
                                                    A=A(type1, type2); \
                                                    rho=rho(type1, type2); \
                                                    C=C(type1, type2)",
                                    "fep" : "lambdaVdw^2*(A/(alpha+r^2)^6 - B/(alpha+r^2)^3)"
                                    })

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

    # Loop over all th custom forcefields defined
    for cft in customForceFields:
        # Check if this forcefield type is used
        isRequired = False
        for iList in interactionsList:
            if iList["ff"] == cft["type"]:
                # isRequired = True
                if (
                    iList["type1"] in listOfAtomTypesPresent
                    and iList["type2"] in listOfAtomTypesPresent
                ):
                    isRequired = True
        if not isRequired:
            continue
        logger.debug("FF: Custom forcefield type = {}".format(cft["type"]))

        # Skip if the entry has no "expression"
        if "expression" not in cft:
            continue
        # Create a new force object
        expression = cft["expression"]
        if setup.FEP and setup.config["fep"]["type"].lower() == "vdw":
            expression = addForcesSwitchFEP(cft["type"], cft["expression"], cft["fep"])
        elif setup.FEP and setup.config["fep"]["type"].lower() == "ff":
            addForcesSwitchFEP2(cft)
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

        # if setup.FEP:
        newForce.addGlobalParameter("lambdaVdw", 1.0)
        newForce.addGlobalParameter("lambdaCoul", 1.0)
        newForce.addGlobalParameter("lambdaFF", 0.0)
        newForce.addGlobalParameter("lambdaFEP", 0.0)

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

        if setup.FEP and setup.config["fep"]["type"].lower() in ["vdw", "ff"]:
            _ = newForce.addPerParticleParameter("fepAtom")
            _ = newForce.addPerParticleParameter("resID")
            for a in listOfAtoms:
                newForce.addParticle([a["itype"], 1, a["resid"]])
        else:
            for a in listOfAtoms:
                newForce.addParticle([a["itype"]])

        # Add the parameters
        for prm in cft["params"]:
            if cft["init"][prm] == 1:
                paramMatrix = np.ones([numTypes, numTypes])
            else:
                paramMatrix = np.zeros([numTypes, numTypes])

            for iList in interactionsList:
                if iList["ff"] == cft["type"]:
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
        newForce.setName(cft["type"] + "_vdwl")
        system.addForce(newForce)


def debugCustomForcefieldInteractions(fparam, listOfAtomTypes, paramMatrix):
    import logging

    logger = logging.getLogger("dynamicEntropy")
    logger.debug("FF: parameter = {}".format(fparam))
    line = "        "
    for xx in range(len(paramMatrix)):
        line += "{:8s} ".format(listOfAtomTypes[xx])
    logger.debug(line)
    for xx in range(len(paramMatrix)):
        line = "{:6s} ".format(listOfAtomTypes[xx])
        for yy in range(len(paramMatrix)):
            x = paramMatrix[xx, yy]
            if x > 0:
                n = int(np.log10(x))
                m = min(6, 6 - n)
                if n > 6:
                    string = "{:<.8g}".format(round(x, m))
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


def addForcesSwitchFEP(name, expr, exprFEP):
    expr = expr.replace(" ", "")
    pattern = re.compile(r"\br\b")
    expr = re.sub(pattern, "d", expr).split(";")

    newExpr = [
        "XX*fep + (1-XX)*ff ; ; fep=" + exprFEP + "; ; ff=" + (";").join(expr) + " ;"
    ]

    # Add this to avoid 1/0
    newExpr.append("d=max(0.1,r)")

    # p = 0 atom is involved in FEP, m is the residue ID
    # XX                  -> 1 means to use FEP force
    # MM = delta(mm1-mm2) -> 1 means same molecules --> no fep
    # PP = delta(pp1*pp2) -> 1 means one of the atoms is perturbed
    # select(x,y,z) = z if x = 0, y otherwise
    # MM == 1 -> returns 0 [NO FEP - same molecule]
    # MM == 0 -> returns PP ->
    #   PP = 0 [FEP] <- 1-delta(0x0) or delta(0x1)
    #   PP = 1 [NO FEP] <- delta(1x1)
    newExpr.append("XX=select(MM,0,PP)")
    newExpr.append("PP=delta(fepAtom1*fepAtom2)")
    newExpr.append("MM=delta(resID1-resID2)")

    return (";").join(newExpr)


def addForcesSwitchFEP2(cft):
    logger = logging.getLogger("dynamicEntropy")

    expr = cft["expression"].split(";")

    if "newexpr" in cft:
        string = cft["newexpr"]["expression"].replace(" ", "")
        newExpr = string.split(";")

        cft["init"].update(cft["newexpr"]["init"])
        cft["params"] += cft["newexpr"]["params"]
    else:
        # Create regex patterns from the list
        # and add "1" to the parameters' names
        string = cft["expression"].replace(" ", "")
        for prm in cft["params"]:
            pattern = re.compile(r"\b" + re.escape(prm) + r"\b")
            string = re.sub(pattern, prm + "1", string)
        newExpr = string.split(";")
        logger.warning("FEP: New FF expression was automatically generated")

        for x in cft["params"]:
            cft["init"][x + "1"] = cft["init"][x]
        cft["params"] += [x + "1" for x in cft["params"]]

    newFF = ["XX*lambdaFF*fep + (1-lambdaFF*XX)*ff"]
    newFF += ["fep=" + newExpr[0]]
    newFF += ["ff=" + expr[0]]
    newFF += newExpr[1:]
    newFF += expr[1:]

    # p = 0 atom is involved in FEP, m is the residue ID
    # XX                  -> 1 means to use FEP force
    # MM = delta(mm1-mm2) -> 1 means same molecules --> no fep
    # PP = delta(pp1*pp2) -> 1 means one of the atoms is perturbed
    # select(x,y,z) = z if x = 0, y otherwise
    # MM == 1 -> returns 0 [NO FEP - same molecule]
    # MM == 0 -> returns PP ->
    #   PP = 0 [FEP] <- 1-delta(0x0) or delta(0x1)
    #   PP = 1 [NO FEP] <- delta(1x1)
    newFF.append("XX=select(MM,0,PP)")
    newFF.append("PP=delta(fepAtom1*fepAtom2)")
    newFF.append("MM=delta(resID1-resID2)")

    cft["expression"] = (";").join(newFF)

    return
