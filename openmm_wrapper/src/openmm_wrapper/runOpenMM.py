#!/usr/bin/env python3
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

########################## import openMM ###############################
import logging
import numpy as np

import openmm as mm
import openmm.app as app
import openmm_wrapper as my


def runOpenMM():
    cmd = my.commandLineParser()
    setup = my.simulationSetup(cmd)
    my.debugDictionary(cmd, name="CMD")

    if cmd["keys"] is not None:
        my.dumpSampleInputFile(cmd, setup)
        quit()
    del cmd["keys"]

    setup.setParameters(cmd)
    # my.debugDictionary(setup.config, name="CONFIG")

    modeller, system = my.createSystem(setup)

    if cmd["minimise"] is not None and cmd["minimise"]:
        if type(cmd["minimise"]) is bool:
            c = {
                "minimise": {"tolerance": 10, "maxIterations": 0, "output": "geopt.pdb"}
            }
            setup.setParameters(c)
        my.minimise(setup, modeller, system)

    elif cmd["energy"] is not None and cmd["energy"]:
        my.singlePointEnergy(setup, modeller, system)

    elif cmd["forwardFluxSampling"] is not None:
        ffs = my.forwardFluxSampling(setup, modeller, system)

    elif cmd["metadynamics"] is not None:
        MTD = my.initialiseMetadynamics(setup, modeller, system)
        simulation = my.initialiseMolecularDynamics(setup, modeller, system)
        my.runMetadynamics(setup, MTD, simulation)

    elif cmd["fep"] is not None:
        my.free_energy_perturbation(setup, modeller, system)

    elif cmd["ti"] is not None:
        my.thermodynamicIntegration(setup, modeller, system)

    elif cmd["md"] is not None:
        if cmd["plumed"] is not None:
            _ = my.createPlumed(cmd["plumed"], system)

        simulation = my.initialiseMolecularDynamics(setup, modeller, system)
        my.runMolecularDynamics(setup, simulation)

    elif cmd["rerun"] is not None:
        my.rerunTrajectory(setup, modeller, system)


if __name__ == "__main__":
    runOpenMM()
