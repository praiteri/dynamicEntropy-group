### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import numpy as np

import openmm as mm
import openmm.app as app
import openmm.unit as unit
import new_openmm_wrapper as my


class simulationSetup(object):
    def __init__(self, cmd):
        self.runID = None
        self.defaultParameters()
        self.setParameters(cmd)

        my.pretty_log({"OpenMM Version": mm.Platform.getOpenMMVersion()}, sep=True)
        my.pretty_log(self.config["basic"])

    ################################################
    def setAtomTypesFromFF(self, ff):
        self.listOfAtomTypesFF = []
        for x in ff._atomTypes:
            self.listOfAtomTypesFF.append(x)

    def getAtomTypesFromFF(self):
        return self.listOfAtomTypesFF

    def createListOfAtoms(self, ff, topology):
        self.listOfAtoms = []
        self.listOfAtomTypesPresent = []

        for res in ff.getMatchingTemplates(topology):
            for af in res.atoms:
                if af.type not in self.listOfAtomTypesPresent:
                    self.listOfAtomTypesPresent.append(af.type)

        bondedToAtom = ff._buildBondedToAtomList(topology)
        for res, r in zip(ff.getMatchingTemplates(topology), topology.residues()):
            _, mapMatches = ff._getResidueTemplateMatches(r, bondedToAtom)

            list1 = [x.index for x in r.atoms()]
            list2 = [a.name for a in res.atoms]
            list3 = [a.type for a in res.atoms]
            for i, j in enumerate(mapMatches):
                self.listOfAtoms.append(
                    {
                        "index": list1[i],
                        "name": list2[j],
                        "type": list3[j],
                        "resid": int(r.index),
                        "itype": self.listOfAtomTypesPresent.index(list3[j]),
                    }
                )
                # "itype" : self.listOfAtomTypesFF.index(af.type)})

    def getListOfAtoms(self):
        return self.listOfAtoms

    def getListOfAtomTypesPresent(self):
        return self.listOfAtomTypesPresent

    def setExclusionsList(self, system, nbonds=3):
        bonds = self.getBondsList()
        listFromBonds = app.forcefield._findExclusions(bonds, nbonds, system.getNumParticles())
        listFromNonBonded = []
        for force in system.getForces():
            if force.getName() == "NonbondedForce":
                for i in range(force.getNumExceptions()):
                    l = force.getExceptionParameters(i)
                    listFromNonBonded.append((l[0], l[1], 0))

        if len(listFromBonds) == len(listFromNonBonded):
            self.listOfExclusions = listFromBonds
        else:
            # This would work unless the FF is OPLS with 1-5 interactions
            # self.listOfExclusions = listFromNonBonded
            exclusions = {(a, b): c for a, b, c in listFromNonBonded}
            exclusions.update({(a, b): c for a, b, c in listFromBonds})
            self.listOfExclusions = [(a, b, c) for (a, b), c in exclusions.items()]
        return

    def getExclusionsList(self):
        return self.listOfExclusions

    def setBondsList(self, topology):
        self.listOfBonds = []
        for b in topology.bonds():
            self.listOfBonds.append((b.atom1.index, b.atom2.index))

    def getBondsList(self):
        return self.listOfBonds

    ################################################
    def setParameters(self, cmd):
        my.replace_entries_recursive(cmd, self.config)

        # Initialise random number generator
        if self.config["input"]["irand"] is None:
            self.rng = np.random.default_rng()
        else:
            self.rng = np.random.default_rng(self.config["input"]["irand"])

        # Define thermal energy
        if not unit.is_quantity(self.config["md"]["temperature"]):
            self.config["md"]["temperature"] *= unit.kelvin
        self.config_kBT = (
            unit.MOLAR_GAS_CONSTANT_R * self.config["md"]["temperature"]
        ).in_units_of(unit.kilojoule_per_mole)
        self.config_beta = 1.0 / self.config_kBT

        # Define configuration for platform
        if self.config["basic"]["Platform"] is None:
            self.config["basic"]["Platform"] = my.checkPlatforms()

        if self.config["basic"]["Platform"] in ["CPU", "Reference"]:
            for x in ["DeviceIndex", "Precision"]:
                self.config["basic"][x] = "Auto"
            self.properties = {}
        elif self.config["basic"]["Platform"] == "OpenCL":
            for x in ["DeviceIndex", "Precision"]:
                self.config["basic"][x] = "Auto"
            self.properties = {}
        elif self.config["basic"]["Platform"] in ["CUDA"]:
            self.properties = {
                x: str(self.config["basic"][x]) for x in ["DeviceIndex", "Precision"]
            }
        elif self.config["basic"]["Platform"] in ["HIP"]:
            self.properties = {
                x: str(self.config["basic"][x]) for x in ["HipDeviceIndex", "Precision"]
            }
        else:
            raise Exception("Platform not recognised")

    ################################################
    def defaultParameters(self):
        self.listOfAtomTypes = ["XX"]

        self.config = {
            "debug": False,
            "basic": {
                "Platform": None,
                "Precision": "mixed",
                "DeviceIndex": 0,
                "HipDeviceIndex": 0,
            },
            "input": {
                "runID": self.runID,
                "irand": None,
                "coordinates": "coord.pdb",
                "forcefield": "forcefield.xml",
                "dry": None,
            },
            "forcefield": {
                "skipScript": False,
                "isAMOEBA": False,
                "isOPLS": False,
                "nonbondedCutoff": 0.9 * unit.nanometer,
                "ewaldErrorTolerance": 1e-5,
                "nonbondedMethod": app.PME,
                "switchDistance": None,
                "customSwitchDistance": None,
                "useDispersionCorrection": False,
                "constraints": None,
                "rigidWater": False,
                "polarization": "mutual",
                "mutualInducedTargetEpsilon": 1e-05,
                "lj14Scale": 0.5,
                "modify": None,
                # "CMMotionRemover": 1000,
            },
            "energy": False,
            "rerun": False,
            "test_forces": False,
            "md": {
                "minimise": False,
                "ensemble": "NVT",
                "thermostat": "LANG",
                "maxDisplacement": None,
                "timestep": 0.001 * unit.picoseconds,
                "numberOfSteps": 10000,
                "simulationTime": None,
                "temperature": 300 * unit.kelvin,
                "thermostatParameter": 1.0 * unit.picoseconds,
                "pressure": 1 * unit.bar,
                "barostat": "ISO",
                "barostatUpdate": 25,
                "reportInterval": 10000,
                "restartFrom": None,
                "CMMotionRemover": {
                    "type": "internal",
                    "reportInterval": 1000,
                    "groups": None,
                },
            },
            "reporters": None,
            "minimise": {"tolerance": 10, "maxIterations": 0, "output": "geopt.pdb"},
            "restraint": None,
            "osmotic": None,
            "catalyst": None,
            "forwardFluxSampling": {
                "lambda": [],
                "nsteps": 1000000,
                "nsample": 1,
                "ntrials": 100000,
                "lambdaID": 0,
                "sigma": 0.1,
                "cv": {
                    "type": "pair",
                    "expression": "0.5*erfc((r-0.31)/0.07)",
                    "set1": [0],
                    "set2": [1, 4, 7, 10, 13, 16],
                    "cutoff": 1.0,
                },
            },
            "fep": {
                "test": False,
                "type": None,
                "select": [],
                "select1": [],
                "select2": [],
                "indices": [],
                "force": [],
                "lambda_sim": None,
                "lambda_eval": [],
                "equilibrationTime": None,
                "equilibrationSteps": None,
                "reportInterval": 1000,
                "minimise": False,
            },
            "metadynamics": {
                "height": 2.5,
                "frequency": 1000,
                "biasFactor": 5,
                "saveFrequency": 1000,
                "biasDir": "./",
                "cv0": None,
                "cv1": None,
                "cv2": None,
            },
            "plumed": None,
            # "ti": {
            #     "force": "tiForce",
            #     "equilibrationSteps": 10000.0,
            #     "output": "ti.{}.out".format(self.runID),
            #     "reportInterval": 1000,
            # },
        }

    ################################################################################
    def dumpParametersFF(self, parms):
        my.pretty_log(
            self.config["input"]["forcefield"],
            title="Forcefield parameters",
            logger="info",
        )

    def dumpParametersMD(self):
        my.pretty_log(self.config["md"], title="Molecular dynamics setup", logger="info")
