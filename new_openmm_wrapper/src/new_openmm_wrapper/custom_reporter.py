### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import numpy as np
import logging
import openmm as mm
import openmm.app as app
import openmm.unit as unit


class forceReporter(object):
    def __init__(self, file, reportInterval):
        self._out = open(file, "w")
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
        system = simulation.context.getSystem()

        # forceGroup = set()
        string = f"{stime:.3f}"
        state = simulation.context.getState(
            getPositions=True, getForces=True, groups=-1
        )
        forces = state.getForces(asNumpy=True)
        forces = forces.value_in_unit(unit.kilojoule_per_mole / unit.nanometer)

        mean = np.mean(forces, axis=0)
        stderr = np.std(forces, axis=0) / np.sqrt(502)
        string += f" {mean[2]:10.3g} {stderr[2]:10.3g}"

        for force in system.getForces():
            name = force.getName()
            forceGroup = force.getForceGroup()

            state = simulation.context.getState(
                getPositions=True, getForces=True, groups=forceGroup
            )
            forces = state.getForces(asNumpy=True)
            forces = forces.value_in_unit(unit.kilojoule_per_mole / unit.nanometer)

            mean = np.mean(forces, axis=0)
            stderr = np.std(forces, axis=0) / np.sqrt(502)
            string += f" {mean[2]:10.3g} {stderr[2]:10.3g}"
            # print(
            #     forceGroup,
            #     name,
            #     mean[2],
            #     stderr[2],
            # )
        # print(string)
        self._out.write(string + "\n")
        self._out.flush()
