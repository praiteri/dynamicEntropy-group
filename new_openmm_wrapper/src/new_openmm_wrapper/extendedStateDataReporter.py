### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import openmm as mm


class ExtendedStateDataReporter(mm.app.StateDataReporter):
    """
    An extension of OpenMM's StateDataReporter_ class, which also outputs pressure if a barostat is present
    """

    # Inspired by https://atomsmm.readthedocs.io/en/latest/_modules/atomsmm/reporters.html#ExtendedStateDataReporter
    # See also https://github.com/openmm/openmm/issues/5129
    def __init__(self, file, reportInterval, **kwargs):
        self._barostat = kwargs.pop("barostat", None)
        super().__init__(file, reportInterval, **kwargs)

    def _constructHeaders(self):
        headers = super()._constructHeaders()
        if self._barostat:
            headers.append("Pressure (bar)")
        return headers

    def _constructReportValues(self, simulation, state):
        values = super()._constructReportValues(simulation, state)
        if self._barostat:
            pressure = self._barostat.computeCurrentPressure(simulation.context)
            values.append(pressure.value_in_unit(mm.unit.bar))
        return values
