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

import openmm as mm
import openmm.app as app
import openmm.unit as unit

import openmm_wrapper as my
import logging


class catchExplosion(object):
    """
    Stores the last nFrames and dumps them if the energy passes a chosen threshold.
    The time separation between the stored frames can be altered using nFreq.
    """

    def __init__(self, settings, context):
        self._reportInterval = settings["nFreq"]
        self._threshold = settings["energyThreshold"]
        self._context = context
        self._nFrames = settings["nFrames"]
        self._pos = [None] * self._nFrames
        self._ncount = 0

    def __del__(self):
        try:
            self._out.close()
        except:
            pass

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, False, True, None)

    def report(self, simulation, state):
        logger = logging.getLogger("dynamicEntropy")
        idx = self._ncount % self._nFrames
        self._pos[idx] = state.getPositions()

        # Energy taken from MD context, this is OK only if the same precision is used in all contexts
        e0 = state.getPotentialEnergy().value_in_unit(unit.kilojoules / unit.mole)

        if e0 > self._threshold:
            self._out = open("boom.pdb", "w")

            for i in range(0, self._nFrames):
                logger.info(
                    "Dumping coordinates for frame {}".format(i - self._nFrames)
                )
                p = self._pos[(idx + 1 + i) % self._nFrames]
                app.PDBFile.writeFile(simulation.topology, p, self._out)

                my.dumpEnergies([self._context], p)

            self._out.flush()
            quit()

        self._ncount += 1
