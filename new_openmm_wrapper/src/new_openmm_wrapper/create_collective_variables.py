### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import copy
import openmm as mm
import openmm.app as app
import new_openmm_wrapper as my


def createCollectiveVariables(cv, topology, system):
    if cv["type"].lower() == "pair":
        collectiveVariablesForce = my.defineNonBondedCV(cv, topology, system)
    elif cv["type"].lower() == "pos":
        collectiveVariablesForce = definePositionCV(cv, topology)
    else:
        raise Exception("Unknown type of CV ({})".format(cv["type"]))
    return collectiveVariablesForce


def definePositionCV(config, topology):
    cvForce = mm.CustomExternalForce(config["expression"])

    if "index1" in config:
        s1 = config["index1"]
    elif "type1" in config:
        s1 = [atom.index for atom in topology.atoms() if atom.name in config["type2"]]
    else:
        raise Exception("Missing index1/type1 in CV definition {}".config)

    assert len(s1) > 0, "No central atoms selected for CV"
    for i in s1:
        cvForce.addParticle(i, [])

    return cvForce


def defineNonBondedCV(config, topology, system):
    # logger = logging.getLogger("dynamicEntropy")

    cvForce = mm.CustomNonbondedForce(config["expression"])
    cvForce.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)

    for _ in range(topology.getNumAtoms()):
        cvForce.addParticle()

    if sum(["index1" in config, "type1" in config]) == 0:
        raise Exception("MTD: index1 or type1 must be used")
    elif sum(["index1" in config, "type1" in config]) == 0:
        raise Exception("MTD: cannot use both index1 and type1 must be used")

    if "index1" in config:
        if type(config["index1"]) is int:
            config["index1"] = str(config["index1"])
        s1 = [eval(x) for x in config["index1"].split(",")]

    elif "type1" in config:
        typeList = config["type1"].split(",")
        s1 = [atom.index for atom in topology.atoms() if atom.name in typeList]
    else:
        raise Exception("Missing index1/type1 in CV definition {}".config)

    if sum(["index2" in config, "type2" in config]) == 0:
        raise Exception("MTD: index1 or type2 must be used")
    elif sum(["index2" in config, "type2" in config]) == 0:
        raise Exception("MTD: cannot use both index2 and type2 must be used")

    if "index2" in config:
        if type(config["index2"]) is int:
            config["index2"] = str(config["index2"])
        s2 = [eval(x) for x in config["index2"].split(",")]

    elif "type2" in config:
        typeList = config["type2"].split(",")
        s2 = [atom.index for atom in topology.atoms() if atom.name in typeList]
    else:
        raise Exception("Missing index2/type2 in CV definition {}".config)

    assert len(s1) > 0, "NBCV - No central atoms selected for CV"
    assert len(s2) > 0, "NBCV - No neighbours atoms selected for CV"

    my.pretty_log(config, title="Defining nonbonded CV", indent=1)
    my.pretty_log({"Number of atoms selected in set 1": len(s1)}, indent=2)
    my.pretty_log(s1, indent=2, logger="debug")
    my.pretty_log({"Number of atoms selected in set 2": len(s2)}, indent=2)
    my.pretty_log(s2, indent=2, logger="debug")

    cvForce.addInteractionGroup(s1, s2)

    cvForce.setCutoffDistance(config["cutoff"])

    return cvForce


def defineCentroidCV():
    """
    CustomCentroidBondForce* force = new CustomCentroidBondForce(2, "0.5*k*distance(g1,g2)^2");
    force->addPerBondParameter("k");
    force->addGroup(particles1);
    force->addGroup(particles2);
    vector<int> bondGroups;
    bondGroups.push_back(0);
    bondGroups.push_back(1);
    vector<double> bondParameters;
    bondParameters.push_back(k);
    force->addBond(bondGroups, bondParameters);
    """


def addWall(system, cvForce, wallExpr, str):
    """
    Add a wall on the CV using an automatically defined function.
    It work only for nonbonded type of CVs.
    """
    wallForce = copy.copy(cvForce)
    wallForce.setEnergyFunction(wallExpr)
    wallForce.setName(str)

    system.addForce(wallForce)
