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
import openmm_wrapper as my

import logging
import numpy as np
from simtk import unit


class FEPReporter(object):
    """
    Class for reporting free energy perturbation data during a simulation.

    Parameters:
    - setup (object): The setup object containing configuration settings.
    - cList (list): List of context objects to compute energies from.
    """

    def __init__(self, setup, cList, **kwargs):
        logger = logging.getLogger("dynamicEntropy")
        self._equilibrationTime = setup.config["fep"][
            "equilibrationTime"
        ].value_in_unit(unit.picosecond)
        logger.debug(
            "FEP: {} = {}".format("equilibration time", self._equilibrationTime)
        )

        self._kbt = setup.config_kBT.value_in_unit(unit.kilojoule_per_mole)
        logger.debug("FEP: {} = {}".format("kBT", self._kbt))

        self._ncount = 0
        self._context = cList  # list of context to compute the energies from
        self._expDE = np.zeros(len(cList))  # sum of exp(-dE/kT) for each context
        self._deltaU0 = []

        self._headers = "{:#>10s} {:^20}".format("# Time", "E0")
        self._format_str = "{:>10.3f} {:-20.10f}"
        for i in range(0, len(cList)):
            self._headers += " {:^20}".format("E" + str(i + 1))
            self._format_str += " {:-20.10f}"
        for i in range(0, len(cList)):
            self._headers += " {:^20}".format("dE" + str(i + 1))
            self._format_str += " {:-20.10f}"
        for i in range(0, len(cList)):
            self._headers += " {:^20}".format("dG" + str(i + 1))
            self._format_str += " {:-20.10f}"
        self._headers += "\n"
        self._format_str += "\n"

        self._out = open(setup.config["fep"]["output"], "w")
        self._out.write(self._headers)
        self._reportInterval = setup.config["fep"]["reportInterval"]

        if "cc" in kwargs:
            self._context_md = kwargs["cc"]

    def __del__(self):
        """
        Closes the output file when the object is deleted.
        """
        self._out.close()

    def describeNextReport(self, simulation):
        """
        Returns a five element tuple.
        The first element is the number of steps until the next report.
        The remaining elements specify whether that report will require
        positions, velocities, forces, and energies respectively.

        Parameters:
        - simulation (object): The simulation object.

        Returns:
        - tuple: A five element tuple specifying the next report details.
        """
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, True, True, None)

    def report(self, simulation, state):
        """
        Reports the free energy perturbation data.

        Parameters:
        - simulation (object): The simulation object.
        - state (object): The state object containing the current simulation state.
        """
        stime = state.getTime().value_in_unit(unit.picosecond)
        if stime < self._equilibrationTime:
            return

        self._ncount += 1
        pos = state.getPositions()
        box = state.getPeriodicBoxVectors()

        # Energy taken from MD context, this is OK only if the same precision is used in all contexts
        e0 = state.getPotentialEnergy().value_in_unit(unit.kilojoules / unit.mole)

        contextEnergy = []
        energyChange = []
        dG = []

        if False:
            my.computeEnergyFromContext(self._context_md)
            for idx, c in enumerate(self._context):
                c = self._context[0]
                c.setPeriodicBoxVectors(box[0], box[1], box[2])
                c.setPositions(pos)
                my.computeEnergyFromContext(c)
            quit()

        for idx, c in enumerate(self._context):
            c.setPeriodicBoxVectors(box[0], box[1], box[2])
            c.setPositions(pos)
            ef = (
                c.getState(getEnergy=True)
                .getPotentialEnergy()
                .value_in_unit(unit.kilojoules / unit.mole)
            )
            contextEnergy.append(ef)
            dU = ef - e0
            energyChange.append(dU)

            # Define a constant shift for the energy difference
            # to avoid over/underflow of the exponential
            # Assume dU at the first evaluation is of the
            # correct order of magnitude
            #
            # \begin{eqnarray}
            # -\beta\Delta G &=& \ln\langle e^{-\beta\Delta U} \rangle \\
            # & =& \ln\Big\langle e^{-\beta\Delta U} \times \displaystyle\frac{e^{\beta\delta}}{e^{\beta\delta}}\Big\rangle \\
            # & =& \ln\frac{\langle e^{-\beta(\Delta U-\delta)} \rangle }{e^{\beta\delta}} \\
            # & =& \ln\langle e^{-\beta(\Delta U-\delta)} \rangle - \beta\delta
            # \end{eqnarray}
            #
            # \begin{equation}
            # \Delta G = \delta - RT\ln\langle e^{-\beta(\Delta U-\delta)} \rangle
            # \end{equation}

            if len(self._deltaU0) < len(self._context):
                self._deltaU0.append(dU)

            dU -= self._deltaU0[idx]
            self._expDE[idx] += np.exp(-(dU) / self._kbt)
            dG.append(
                self._deltaU0[idx] - self._kbt * np.log(self._expDE[idx] / self._ncount)
            )

        self._out.write(
            self._format_str.format(stime, e0, *contextEnergy, *energyChange, *dG)
        )
        self._out.flush()
