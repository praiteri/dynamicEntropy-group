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

import argparse
import pathlib
import yaml
import json

import sys
import openmm_wrapper as my
import pprint

pp = pprint.PrettyPrinter(indent=4)


class ParseDict(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        d = getattr(namespace, self.dest) or {}
        # to cope with commands where there are a spaces around the "=" sign
        xx = []
        for i, v in enumerate(values):
            if len(v) > 1 and "=" in v:
                xx.append(v)
            else:
                if v == "=":
                    xx.append("".join(values[i - 1 : i + 2]))
        # values = [ "".join(values[i-1:i+2]) for i,v in enumerate(values) if v == "=" ]
        if xx:
            for item in xx:
                split_items = item.split("=", 1)
                key = split_items[0].strip()  # we remove blanks around keys, as is logical
                value = split_items[1]
                d[key] = value
        setattr(namespace, self.dest, d)


def commandLineParser():
    """
    Parse the command line arguments and put them in a dictionary.
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--debug",
        dest="debug",
        type=bool,
        required=False,
        action=argparse.BooleanOptionalAction,
        help="""
Provide verbose debugging output.
                        """,
    )

    parser.add_argument(
        "--log",
        dest="log",
        type=str,
        required=False,
        action="append",
        help="""
Provide verbose debugging output.
                        """,
    )

    parser.add_argument(
        "--b",
        "--basic",
        dest="basic",
        type=str,
        required=False,
        metavar="KEY=VALUE",
        nargs="*",
        action=ParseDict,
        help="""
This command can be used to set low levwl hardware options.
Available KEY options and default values:
    "Platform" : "CUDA"
    "Precision" : "mixed"
    "DeviceIndex" : "0"
                        """,
    )

    parser.add_argument(
        "--i",
        "--input",
        dest="input",
        type=str,
        required=False,
        metavar="KEY=VALUE",
        nargs="*",
        action=ParseDict,
        help="""
This command can be used to set the simulation's coordinates and forcefield input files.
Available KEY options and default values:
    "runID"                      : automatic counter based on the number of output.*.out files
    "irand"                      : None
    "coordinates"                : "coord.pdb"
    "forcefield"                 : "newFF.xml"
    "ensemble"                   : "NVT",
                        """,
    )

    parser.add_argument(
        "--ff",
        "--forcefield",
        dest="forcefield",
        type=str,
        required=False,
        metavar="KEY=VALUE",
        nargs="*",
        action=ParseDict,
        help="""
This command can be used to set the forcefield global parameters.
Available KEY options and default values:
    "isAMOEBA"                   : False,
    "isOPLS"                     : False,
    "nonbondedCutoff"            : 0.9 * unit.nanometer,
    "ewaldErrorTolerance"        : 1e-5,
    "nonbondedMethod"            : app.PME,
    "switchDistance"             : None,
    "customSwitchDistance"       : None,
    "useDispersionCorrection"    : False,
    "constraints"                : None,
    "rigidWater"                 : False,
    "polarization"               : "mutual",
    "mutualInducedTargetEpsilon" : 1e-05,
    "lj14Scale"                  : 0.5,
                        """,
    )

    parser.add_argument(
        "--e",
        "--energy",
        dest="energy",
        type=bool,
        required=False,
        action=argparse.BooleanOptionalAction,
        help="""
This command does a single point calculation and write the energy breakdown.
                        """,
    )

    parser.add_argument(
        "--minimise",
        dest="minimise",
        type=str,
        metavar="KEY=VALUE",
        nargs="*",
        action=ParseDict,
        help="""
This command performs an energy minimisation.
Available KEY options and default values:
    "tolerance"                  : 10
    "maxIterations"              : 0
    "output"                     : "geopt.pdb"                        
                        """,
    )

    parser.add_argument(
        "--md",
        "--molecularDynamics",
        dest="md",
        type=str,
        metavar="KEY=VALUE",
        nargs="*",
        action=ParseDict,
        help="""
This command can be used to set the parameters for the Molecular Dynamics part of a simulation.
Available KEY options and default values:
    "ensemble"                   : "NVT"
    "thermostat"                 : "LANG"
    "timestep"                   : 0.001 (unit.picoseconds)
    "numberOfSteps"              : 10000
    "simulationTime"             : -1
    "temperature"                : 300 (unit.kelvin)
    "thermostatParameter"        : 1.0 (unit.picoseconds or 1/unit.picoseconds)
    "pressure"                   : 1 (unit.bar)
    "barostat"                   : "ISO"
    "barostatUpdate"             : 25
    "reportInterval"             : 1000
    "restartFrom"                : ""
                        """,
    )

    parser.add_argument(
        "--r",
        "--restraint",
        dest="restraint",
        type=str,
        required=False,
        metavar="KEY=VALUE",
        nargs="*",
        action="append",
        help="""
This command adds positional restratins to the atoms.
    Sample restraint commands:
        --restraint '{ "r1": {"file":"rest_bot.pdb","par":["k"], "val":[10000]}}' 
        --restraint '{ "r2" : {"par":["k","d0"],"val":[10000,5],"global":"d0","species":["Na","Cl"],"fullexp":"k*max(0,y-d0)^2" } }' 
        --restraint '{ "r3" : {"par":["k","d0"],"val":[0,5],"global":"d0","species":["Na","Cl"],"fullexp":"k*max(0,y-d0)^2" } }'                        
    The restraint command can be repeated on the command line.
    In the YAML file the restraints should be defined as a dictionary of dictionaries.
        """,
    )
    parser.add_argument(
        "--osm",
        "--osmotic",
        dest="osmotic",
        type=str,
        required=False,
        metavar="KEY=VALUE",
        nargs="*",
        action="append",
        help="""
                        """,
    )

    parser.add_argument(
        "--ffs",
        "--forwardFlux",
        dest="forwardFluxSampling",
        type=str,
        required=False,
        metavar="KEY=VALUE",
        nargs="*",
        action=ParseDict,
        help="""
This command can be used to run Forward Flux Sampling simulations.
Available KEY options and default values:
    "lambda"                     : []
    "nsteps"                     : 1000000
    "nsample"                    : 10
    "ntrials"                    : 100000
    "lambdaID"                   : 0
    "cv"                         : { "type"       : "nb"
                                   , "expression" : "0.5*erfc((r-0.31)/0.07)"
                                   , "set1"       : [0]
                                   , "set2"       : [1,4,7,10,13,16]
                                   , "cutoff"     : 1.0
                        """,
    )

    parser.add_argument(
        "--mtd",
        "--metadynamics",
        dest="metadynamics",
        type=str,
        required=False,
        metavar="KEY=VALUE",
        nargs="*",
        action=ParseDict,
        help="""
                        """,
    )

    parser.add_argument(
        "--plumed",
        dest="plumed",
        type=str,
        required=False,
        metavar="KEY=VALUE",
        nargs="*",
        help="""
                        """,
    )

    parser.add_argument(
        "--fep",
        "--free_energy_perturbation",
        dest="fep",
        type=str,
        required=False,
        metavar="KEY=VALUE",
        nargs="*",
        action=ParseDict,
        help="""
                        """,
    )

    parser.add_argument(
        "--ti",
        "--thermodynamicIntegration",
        dest="ti",
        type=str,
        metavar="KEY=VALUE",
        nargs="*",
        action=ParseDict,
        help="""
This command performs thermodynamic integration calculation.
    Examples of valid options are:
                        """,
    )

    parser.add_argument(
        "--rerun",
        "--rerun",
        dest="rerun",
        type=str,
        required=False,
        metavar="KEY=VALUE",
        nargs="*",
        action=ParseDict,
        help="""
                        """,
    )

    #########################
    #     parser.add_argument("--so", "--screenOutput", dest="screenOutput", type=str, required=False,
    #                         metavar="KEY=VALUE", nargs="*",
    #                         action=ParseDict,
    #                         help="""
    # This command can be used to modify the screen output.
    # Available KEY options and default values:
    #     "file"                       : stdout
    #     "reportInterval"             : 1000
    #     "separator"                  : " "
    #     "step"                       : False
    #     "time"                       : False
    #     "potentialEnergy"            : True
    #     "kineticEnergy"              : False
    #     "totalEnergy"                : False
    #     "temperature"                : True
    #     "volume"                     : True
    #     "density"                    : False
    #     "progress"                   : True
    #     "remainingTime"              : True
    #     "speed"                      : True
    #     "elapsedTime"                : True
    #                         """)

    #     parser.add_argument("--lo", "--logOutput", dest="logOutput", type=str, required=False,
    #                         metavar="KEY=VALUE", nargs="*",
    #                         action=ParseDict,
    #                         help="""
    # This command can be used to modify the screen output.
    # Available KEY options and default values:
    #     "file"                       : "output.{runID}.out"
    #     "reportInterval"             : 1000
    #     "separator"                  : ""
    #     "step"                       : False
    #     "time"                       : True
    #     "potentialEnergy"            : True
    #     "kineticEnergy"              : False
    #     "totalEnergy"                : False
    #     "temperature"                : True
    #     "volume"                     : True
    #     "density"                    : True
    #     "progress"                   : False
    #     "remainingTime"              : False
    #     "speed"                      : True
    #     "elapsedTime"                : True
    #                         """)

    #     # parser.add_argument("--to", "--trajectoryOutput", dest="trajectoryOutput", type=str, required=False,
    #     #                     metavar="KEY=VALUE", nargs="*",
    #     #                     action=ParseDict,
    #     #                     help="""
    #     # "file"                       : "trajectory.{runID}.dcd"
    #     # "reportInterval"             : 1000
    #     # "enforcePeriodicBox"         : True
    #     #                     """)

    #     parser.add_argument("--ro", "--restartOutput", dest="restartOutput", type=str, required=False,
    #                         metavar="KEY=VALUE", nargs="*",
    #                         action=ParseDict,
    #                         help="""
    #     "file"                       : "restart.{runID}.xml"
    #     "checkpointInterval"         : 1000000
    #     "checkpointFile"             : "state.chk"

    #                         """)

    parser.add_argument(
        "-keys",
        "--keys",
        dest="keys",
        type=str,
        required=False,
        metavar="KEY=VALUE",
        nargs="*",
        action="store",
        help="""
    Writes a sample YAML input file
                        """,
    )

    # Parse the command line
    args, unknown = parser.parse_known_args()
    args = vars(args)

    for u in unknown:
        # Read YAML file
        if pathlib.Path(u).suffix in [".yaml", ".yml"]:
            with open(u, "r") as f:
                cmdDict = yaml.safe_load(f)

        # Read JSON file
        elif pathlib.Path(u).suffix in [".json"]:
            with open(u, "r") as f:
                cmdDict = json.load(f)

        # Throw an exception for unknown command
        else:
            raise Exception("Uknown command ({})".format(u))

        # Merge the command line arguments with those in the file
        for k, v in cmdDict.items():
            if k not in args:
                # raise ValueError("command {} not recognised".format(k))
                args[k] = v
            elif args[k] is None:
                args[k] = v
            else:
                try:
                    args[k] = {**v, **args[k]}
                except:
                    pass

    return args
