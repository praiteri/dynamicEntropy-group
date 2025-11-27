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

import numpy as np
import logging
import openmm as mm
import openmm.app as app
import openmm.unit as unit


class forceReporter(object):
    def __init__(self, file, reportInterval):
        self._out = open(file, "w")
        self._out2 = open("pos.dat", "w")
        self._reportInterval = reportInterval

    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        """
        # Returns a five element tuple.
        # The first element is the number of steps until the next report.
        # The remaining elements specify whether that report will require
        # positions, velocities, forces, and energies respectively
        """
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, True, False, None)

    def report(self, simulation, state):
        stime = state.getTime().value_in_unit(unit.picosecond)

        forceGroup = set()
        system = simulation.context.getSystem()
        for force in system.getForces():
            if force.getName() != "CustomCVForce":
                forceGroup.add(force.getForceGroup())

        state = simulation.context.getState(getPositions=True, getForces=True, groups=forceGroup)
        forces = state.getForces(asNumpy=True)
        forces = forces.value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
        positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)

        self._out.write("#Simulation time = {}\n".format(stime))
        for i in range(0, 4):
            self._out.write("%g %g %g " % (forces[i][0], forces[i][1], forces[i][2]))
        self._out.write("\n")
        self._out.flush()

        self._out2.write("#Simulation time = {}\n".format(stime))
        for i in range(0, 4):
            self._out2.write("%g %g %g " % (positions[i][0], positions[i][1], positions[i][2]))
        self._out2.write("\n")
        self._out2.flush()
