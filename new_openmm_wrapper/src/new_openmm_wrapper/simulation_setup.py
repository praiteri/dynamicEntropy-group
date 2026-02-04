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
        if cmd is None:
            return

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
            "minimise": {"tolerance": 10, "maxIterations": 0, "output": "geopt.pdb"},
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
            "restraint": None,
            "osmotic": None,
            "catalyst": None,
            "fep": {
                "path": ".",
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
            # "rerun": False,
            # "test_forces": False,
            # "ti": {
            #     "force": "tiForce",
            #     "equilibrationSteps": 10000.0,
            #     "output": "ti.{}.out".format(self.runID),
            #     "reportInterval": 1000,
            # },
            # "forwardFluxSampling": {
            #     "lambda": [],
            #     "nsteps": 1000000,
            #     "nsample": 1,
            #     "ntrials": 100000,
            #     "lambdaID": 0,
            #     "sigma": 0.1,
            #     "cv": {
            #         "type": "pair",
            #         "expression": "0.5*erfc((r-0.31)/0.07)",
            #         "set1": [0],
            #         "set2": [1, 4, 7, 10, 13, 16],
            #         "cutoff": 1.0,
            #     },
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


def write_sample_input(args):
    import io
    import yaml
    import sys
    import numpy as np
    import openmm.unit as unit
    import openmm.app as app

    def fix_entries(my_dict):
        if not isinstance(my_dict, dict):
            return my_dict

        for k, v in my_dict.items():
            if isinstance(v, dict):
                fix_entries(v)
            else:
                if isinstance(v, type(app.forcefield.PME)):
                    my_dict[k] = f"{v}"
                elif unit.is_quantity(v):
                    my_dict[k] = v._value
                elif isinstance(v, io.TextIOWrapper):
                    my_dict[k] = "stdout"
                else:
                    my_dict[k] = v

        return my_dict

    setup = simulationSetup(None)
    output_dict = fix_entries(setup.config)

    minimal_output = True
    if len(args) == 1 and "," in args[0]:
        args = args[0].split(",")

    if len(args) == 0 or (len(args) == 1 and args[0] == "all"):
        list_of_fields = {k: 1 for k in setup.config}

    else:
        list_of_fields = {k: 0 for k in setup.config}

        for item in args:
            if item in list_of_fields:
                list_of_fields[item] = 1
            else:
                if item == "full":
                    minimal_output = False
                elif item == "all":
                    list_of_fields = {k: 1 for k in setup.config}
                else:
                    my.pretty_log(f"Unknown input section: {item}", logger="error")
                    sys.exit(1)

    if not minimal_output:
        if sum(list_of_fields.values()) == 0:
            list_of_fields = {k: 1 for k in setup.config}
        else:
            for x in [
                "debug",
                "basic",
                "input",
                "forcefield",
            ]:
                list_of_fields[x] = 1

    log_str = [k for k, v in list_of_fields.items() if v == 1]
    my.pretty_log(
        log_str,
        title="Example input file for:",
        sep=True,
        align_width=0,
        ncols=4,
    )

    to_remove = []
    for k in output_dict.keys():
        if list_of_fields[k] == 0:
            to_remove.append(k)
    for k in to_remove:
        del output_dict[k]

    if list_of_fields["minimise"] == 1:
        if not minimal_output:
            output_dict["minimise"] = fix_entries(setup.config["minimise"])
        else:
            output_dict["minimise"] = True

    if list_of_fields["md"] == 1:
        if not minimal_output:
            list_of_fields["reporters"] = 1

    if list_of_fields["reporters"] == 1:
        if not minimal_output:
            output_dict["reporters"] = fix_entries(my.simulationReporters().reporters)
        else:
            output_dict["reporters"] = {
                "screen": True,
                "log": True,
                "dcd": False,
                "xtc": {"file": "traj.xtc", "reportInterval": 2000},
                "restart": {
                    "xml": {"file": "restart.xml", "reportInterval": 100000},
                    "chk": {"file": "restart.chk", "reportInterval": 100000},
                },
            }

    if "restraint" in list(output_dict.keys()):
        output_dict["restraint"] = example_restraint_cmd()

    # else:
    #     list_of_fields = set()
    #     for item in args:
    #         if "," in item:
    #             list_of_fields += item.split(",")
    #         else:
    #             list_of_fields.append(item)

    # for item in list_of_fields:
    #     if item == "md":
    #         list_of_fields +=

    my.pretty_log(sep=True)
    yaml.dump(output_dict, sys.stdout, default_flow_style=False, sort_keys=False, indent=4)
    sys.exit(0)


def example_restraint_cmd():
    return {
        "wall_b": {
            "fullexp": "0.5*k0*(max(0,d-d0))^2; d=abs(z-z0)",
            "par": "k0,z0,d0",
            "select": {"species": "O1"},
            "val": "10000,9.0,2.0",
        },
        "restraintForce": {
            "fullexp": "0.5*k*d^2; d = periodicdistance(x,y,z,x0,y0,z0)",
            "par": "k,x0,y0,z0",
            "val": 10000,
            "file": "rest.pdb",
        },
    }
