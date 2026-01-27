### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import sys
import new_openmm_wrapper as my
import openmm.app as app


class simulationReporters:
    """Default reporters configuration."""

    def __init__(
        self,
        simulation,
        runID=0,
        reportInterval=1000,
        numberOfSteps=None,
        configReporters=None,
    ):
        self.runID = runID
        self.reportInterval = reportInterval

        self.reporters = {
            "screen": {
                "file": sys.stdout,
                "reportInterval": 0,
                "separator": "\t",
                "step": False,
                "time": True,
                "potentialEnergy": True,
                "kineticEnergy": False,
                "totalEnergy": False,
                "temperature": True,
                "volume": True,
                "density": False,
                "progress": False,
                "remainingTime": False,
                "speed": True,
                "elapsedTime": True,
                "barostat": None,
                "totalSteps": None,
            },
            "log": {
                "file": "output.{}.out".format(self.runID),
                "reportInterval": self.reportInterval,
                "separator": ",",
                "step": False,
                "time": True,
                "potentialEnergy": True,
                "kineticEnergy": True,
                "totalEnergy": False,
                "temperature": True,
                "volume": True,
                "density": True,
                "progress": False,
                "remainingTime": False,
                "speed": True,
                "elapsedTime": True,
                "barostat": None,
            },
            "dcd": {
                "file": "trajectory.{}.dcd".format(self.runID),
                "reportInterval": self.reportInterval,
                "enforcePeriodicBox": True,
            },
            "xtc": {
                "file": "trajectory.{}.xtc".format(self.runID),
                "reportInterval": self.reportInterval,
                "enforcePeriodicBox": True,
            },
            "restart": {
                "file": "restart.{}.xml".format(self.runID),
                "reportInterval": 1000000,
                "checkpoint": "state.chk",
            },
        }

        if numberOfSteps is not None:
            self.reporters["screen"]["reportInterval"] = max(1, numberOfSteps // 100)
            self.reporters["screen"]["totalSteps"] = numberOfSteps
            self.reporters["screen"]["progress"] = True
            self.reporters["screen"]["remainingTime"] = True

        self.active_reporters = ["screen", "log"]

        if configReporters is not None:
            for rep, val in configReporters.items():
                if not val and rep in self.active_reporters:
                    self.active_reporters.remove(rep)
                elif val and rep not in self.active_reporters:
                    self.active_reporters.append(rep)

            for k, v in configReporters.items():
                if k in self.reporters and isinstance(v, dict):
                    self.reporters[k].update(v)

        if "xtc" in self.active_reporters and "dcd" in self.active_reporters:
            my.pretty_log(
                "Both XTC and DCD reporters are active",
                logger="WARNING",
            )

        for rep in self.active_reporters:
            my.pretty_log(f" Creating {rep.upper()} reporter", sep=True)

            if rep == "screen":
                self.addScreenOutput(simulation, self.reporters[rep])
            elif rep == "log":
                self.addFileOutput(simulation, self.reporters[rep])
            elif rep == "dcd":
                self.addDCDOutput(simulation, self.reporters[rep])
            elif rep == "xtc":
                self.addXTCOutput(simulation, self.reporters[rep])

        return

    def addScreenOutput(self, simulation, config):
        my.pretty_log(
            config,
            title="Screen reporter configuration:",
            indent=1,
            logger="debug",
        )
        simulation.reporters.append(my.ExtendedStateDataReporter(**config))

    def addFileOutput(self, simulation, config):
        my.pretty_log(
            config,
            title="Logfile reporter configuration:",
            indent=1,
            logger="debug",
        )
        simulation.reporters.append(my.ExtendedStateDataReporter(**config))

    def addDCDOutput(self, simulation, config):
        my.pretty_log(
            config,
            title="DCD reporter configuration:",
            indent=1,
            logger="debug",
        )
        simulation.reporters.append(app.DCDReporter(**config))

    def addXTCOutput(self, simulation, config):
        my.pretty_log(
            config,
            title="XTC reporter configuration:",
            indent=1,
            logger="debug",
        )
        simulation.reporters.append(app.XTCReporter(**config))
