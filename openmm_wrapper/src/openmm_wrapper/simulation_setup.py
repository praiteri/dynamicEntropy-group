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
from sys import stdout

import openmm as mm
import openmm.app as app
import openmm.unit as unit
import openmm_wrapper as my

import subprocess
import numpy as np
import pprint as pp


class simulationSetup(object):
    def __init__(self, cmd):
        self.debug = False

        try:
            self.runID = str(cmd["input"]["runID"])
        except:
            runID = subprocess.getoutput("ls output.*.out 2> /dev/null | wc -l")
            self.runID = runID.strip()

        # Set up logger
        logger_level = "INFO"
        if cmd["debug"] is not None:
            self.debug = True
            logger_level = "DEBUG"
        del cmd["debug"]

        logger_file = None
        if cmd["log"] is not None:
            logger_file = cmd["log"][0]
            # Remove screen output from MD reporters
            if cmd["md"] is not None:
                try:
                    cmd["md"]["screenOutput"]["file"] = None
                except:
                    cmd["md"]["screenOutput"] = {"file": None}
        del cmd["log"]

        logger = my.setup_logger(
            name="dynamicEntropy", level=logger_level, filename=logger_file
        )
        logger.critical("#------------------------------------------#")

        # Set up default configuration parameters
        self.defaultParameters()

        self.config["basic"]["Platform"] = my.checkPlatforms()

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
        listFromBonds = app.forcefield._findExclusions(
            bonds, nbonds, system.getNumParticles()
        )
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
        # Scree output
        self.dumpBasicInfo()

        self.FEP = False
        if "fep" in cmd:
            if cmd["fep"]:
                self.FEP = True

    ################################################
    def defaultParameters(self):
        self.listOfAtomTypes = ["XX"]

        self.config = {
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
                "trajectoryOutput": {
                    "file": "trajectory.{}.xtc".format(self.runID),
                    "reportInterval": 0,
                    "enforcePeriodicBox": True,
                    # "removeWater"                : False,
                    # "removeVirtualSites"         : True,
                },
                "screenOutput": {
                    "file": "stdout",
                    "reportInterval": 0,
                    "separator": "\t",
                    "step": False,
                    "time": True,
                    "potentialEnergy": True,
                    "kineticEnergy": False,
                    "totalEnergy": False,
                    "temperature": True,
                    "volume": True,
                    "density": False,
                    "progress": True,
                    "remainingTime": True,
                    "speed": True,
                    "elapsedTime": True,
                    "barostat": None,
                },
                "logOutput": {
                    "file": "output.{}.out".format(self.runID),
                    "reportInterval": 0,
                    "separator": ",",
                    "step": False,
                    "time": True,
                    "potentialEnergy": True,
                    "kineticEnergy": True,
                    "totalEnergy": False,
                    "temperature": True,
                    "volume": True,
                    "density": True,
                    "progress": False,
                    "remainingTime": False,
                    "speed": True,
                    "elapsedTime": True,
                    "barostat": None,
                },
                "restartOutput": {
                    "file": "restart.{}.xml".format(self.runID),
                    "checkpointInterval": 1000000,
                    "checkpointFile": "state.chk",
                },
            },
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
                "indices": [],
                "species": [],
                "residue": [],
                "force": [],
                "lambda": [],
                "equilibrationTime": 10 * unit.picoseconds,
                "output": "fep.{}.out".format(self.runID),
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
            "ti": {
                "force": "tiForce",
                "equilibrationSteps": 10000.0,
                "output": "ti.{}.out".format(self.runID),
                "reportInterval": 1000,
            },
        }

    ################################################################################
    def dumpParametersFF(self, parms):
        logger = logging.getLogger("dynamicEntropy")
        logger.critical("Opening forcefield ...")
        logger.info("  {:40s} = {}".format("file", self.config["input"]["forcefield"]))
        for x in parms.items():
            logger.info("  {:40s} = {}".format(*x))

    def dumpParametersMD(self):
        logger = logging.getLogger("dynamicEntropy")
        logger.critical("Molecular dynamics setup ...")
        level = logging.getLevelName(logger.getEffectiveLevel())
        for x in self.config["md"].items():
            if isinstance(x[1], dict):
                if level == "DEBUG":
                    logger.info("  {:40s}".format(x[0]))
                    for k, v in x[1].items():
                        logger.debug("   {:<20s} : {}".format(k, v))
                else:
                    if x[1].get("file", None) is not None:
                        logger.info(
                            "  {:40s} = {}".format(x[0] + ".file", x[1]["file"])
                        )
                    else:
                        for k, v in x[1].items():
                            logger.info("  {:40s} = {}".format(k, v))
            else:
                logger.info("  {:40s} = {}".format(*x))
        return

    def dumpBasicInfo(self):
        logger = logging.getLogger("dynamicEntropy")
        logger.critical("openMM parameters ...")
        logger.info("  {:40s} = {}".format("Version", mm.Platform.getOpenMMVersion()))
        for x in self.config["basic"].items():
            if x[0] == "Platform" and x[1] not in ["CUDA", "HIP"]:
                logger.warning("  {:40s} = {}".format(*x))
            else:
                logger.info("  {:40s} = {}".format(*x))
        logger.info("  {:40s} = {}".format("runID", self.runID))
        return


def dumpInfo(string, config):
    logger = logging.getLogger("dynamicEntropy")
    logger.info("  {:40s}".format(string + ":"))
    for x in config.items():
        logger.info("    {:38s} = {}".format(*x))
    return


def cleanUnitsForYaml(dictionary):
    for key, value in dictionary.items():
        if isinstance(value, dict):
            cleanUnitsForYaml(value)
        else:
            if isinstance(value, unit.quantity.Quantity):
                dictionary[key] = value._value
            else:
                dictionary[key] = "{}".format(value)


def addSampleCommands(cmd):
    cmd["restraint"] = {
        "fileRestraint": {
            "file": "rest_bot.pdb",
            "par": "k",
            "val": 10000,
        },
        "wall": {
            "exp": "0.5*k*min(0,z-z0)^2",
            "par": "k,z0",
            "val": "100,45",
            "species": "Ar",
        },
    }

    cmd["osmotic"] = {
        "reportInterval": 1000,
        "computeInterval": 1000,
        "sampleLength": 10,
        "force": "osmoticWall",
        "gcmd": {"pext": 1.0, "tau": 1.0},
        "geometry": {
            "sphere": {
                "centre": "45,45,45",
                "radius": 22.5,
                "kappa": 100.0,
                "species": "Ar",
            },
            "cylinder": {
                "axis": "y",
                "centre": "45,45",
                "kappa": 100,
                "species": "Ar",
                "radius": 20,
            },
            "slab": {
                "axis": "z",
                "width": 22.5,
                "centre": 45,
                "kappa": 100,
                "species": "Ar",
            },
            "plane": {"axis": "-x", "pos": 7, "kappa": 1000.0, "species": "Na,Cl"},
        },
    }

    return


def dumpSampleInputFile(cmd, setup):
    import yaml

    logger = logging.getLogger("dynamicEntropy")

    cleanUnitsForYaml(setup.config)
    addSampleCommands(setup.config)

    requiredFields = [
        "basic",
        "input",
        "forcefield",
        "energy",
        # "rerun",
        "md",
        "minimise",
        # "restraint",
        # "osmotic",
        # "catalyst",
        # "forwardFluxSampling",
        # "fep",
        # "metadynamics",
        # "plumed",
        # "ti",
    ]
    for x in cmd["keys"]:
        if x not in setup.config:
            options = [y for y in setup.config]
            raise Exception(
                "Available options for --keys are: {}".format(", ".join(options))
            )
        requiredFields.append(x)

    dumpDict = {}
    for x in setup.config:
        if x in requiredFields:
            if setup.config[x] is not None:
                dumpDict[x] = setup.config[x]
            else:
                dumpDict[x] = addSampleCommands(x)

    fout = "de2.yaml"
    logger.info(f"Writing sample input file for '{cmd['keys']}' to '{fout}'")
    with open(fout, "w") as f:
        yaml.dump(dumpDict, f, default_flow_style=False, indent=4)
