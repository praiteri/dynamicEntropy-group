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


class tiReporter(object):
    def __init__(self, setup):
        outFile = setup.config["ti"]["output"]
        self._out = open(outFile, "w")
        self._reportInterval = setup.config["ti"]["reportInterval"]
        self._name = setup.config["ti"]["force"]

    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, True, None)
        # A five element tuple.
        # The first element is the number of steps until the next report.
        # The remaining elements specify whether that report will require
        # positions, velocities, forces, and energies respectively

    def report(self, simulation, state):
        forceGroup = set()
        system = simulation.context.getSystem()
        for force in system.getForces():
            if force.getName() == self._name:
                forceGroup.add(force.getForceGroup())
                break

        # extract forces
        l = simulation.context.getParameter("lambdaTI")
        if l > 0:
            erg = simulation.context.getState(
                getEnergy=True, groups=forceGroup
            ).getPotentialEnergy()
            erg = erg.value_in_unit(unit.kilojoule_per_mole) / l
        else:
            erg = 0

        self._out.write("%g %g %g\n" % (l, erg, erg * l))
        self._out.flush()


def thermodynamicIntegration(setup, modeller, system):
    logger = logging.getLogger("dynamicEntropy")

    # Create simulation object
    simulation = my.initialiseMolecularDynamics(setup, modeller, system)

    # Thermodynamic intergration reporter
    simulation.reporters.append(my.tiReporter(setup))

    # Set intial value for the lambda parameter
    c = simulation.context
    c.setParameter("lambdaTI", 0)

    # Equilibrate the system
    logger.critical("#--- System equilibration l=0")
    simulation.step(setup.config["ti"]["equilibrationSteps"])

    # Number of MD stepd for each TI run
    n = int(setup.config["md"]["numberOfSteps"])

    # Run thermodynamic integration l=0 -> l=1
    logger.critical("#--- Forward TI (l=0 -> l=1)")
    for i in range(1, n + 1):
        c.setParameter("lambdaTI", i / n)
        simulation.step(1)

    # Run thermodynamic integration
    logger.critical("#--- System equilibration l=1")
    simulation.step(setup.config["ti"]["equilibrationSteps"])

    # Run thermodynamic integration l=1 -> l=0
    logger.critical("#--- Backward TI (l=1 -> l=0)")
    for i in range(n, 0, -1):
        c.setParameter("lambdaTI", i / n)
        simulation.step(1)
