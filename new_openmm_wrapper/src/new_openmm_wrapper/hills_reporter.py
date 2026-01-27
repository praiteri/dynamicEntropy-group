### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import numpy as np
import openmm.app as app
import openmm.unit as unit


class HILLSReporter(object):

    def __init__(self, file, reportInterval, meta: app.Metadynamics):
        self._out = open(file, "a")
        self._reportInterval = reportInterval
        self.meta = meta
        if self.meta.frequency != reportInterval:
            raise ValueError("HILL frequency must be same as Meta frequency")
        self.gaussian_widths = [v.biasWidth for v in self.meta.variables]
        self.hills_fmt = self.init_header()

    def init_header(self) -> str:
        n_cvs = len(self.meta.variables)
        str_cv = str_s = str_f1 = str_f2 = ""
        for i in range(n_cvs):
            str_cv += " cv" + str(i)
            str_s += " sigma_cv" + str(i)
            str_f1 += " {1[" + str(i) + "]:<16.8f}"
            str_f2 += " {" + str(i + 2) + ":<16.8f}"
        fmt_str = (
            "{0:<16.8f}"
            + str_f1
            + str_f2
            + " {"
            + str(n_cvs + 2)
            + ":<16.8f} {"
            + str(n_cvs + 3)
            + ":<16}\n"
        )
        self._out.write(f"#! FIELDS time{str_cv} {str_s} height biasf\n")
        return fmt_str

    def describeNextReport(self, simulation):
        """
        # Returns a five element tuple.
        # The first element is the number of steps until the next report.
        # The remaining elements specify whether that report will require
        # positions, velocities, forces, and energies respectively
        """
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, True, False, False, None)

    def report(self, simulation, state):
        self.meta._force.updateParametersInContext(simulation.context)

        position = self.meta._force.getCollectiveVariableValues(simulation.context)
        energy = simulation.context.getState(
            getEnergy=True, groups={self.meta._force.getForceGroup()}
        ).getPotentialEnergy()
        height = self.meta.height * np.exp(
            -energy / (unit.MOLAR_GAS_CONSTANT_R * self.meta._deltaT)
        )
        time = simulation.context.getTime().value_in_unit(unit.picosecond)
        scaled_height = (self.meta.biasFactor / (self.meta.biasFactor - 1)) * height
        self._out.write(
            self.hills_fmt.format(
                time,
                position,
                *self.gaussian_widths,
                scaled_height._value,
                self.meta.biasFactor,
            )
        )
        self._out.flush()

    def __del__(self):
        if hasattr(self, "_out"):
            self._out.close()
