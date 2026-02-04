### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import sys
import new_openmm_wrapper as my
import openmm.app as app
import glob


class simulationReporters:
    """Default reporters configuration."""

    def __init__(
        self,
        simulation=None,
        runID=None,
        reportInterval=1000,
        numberOfSteps=None,
        configReporters=None,
    ):

        self.runID = runID
        if runID is None:
            self.runID_str = "."
        else:
            if isinstance(runID, str) and runID.lower() == "auto":
                runID = len(glob.glob("output.*.out"))
            self.runID_str = f".{runID}."

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
                "file": "output{}out".format(self.runID_str),
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
                "file": "trajectory{}dcd".format(self.runID_str),
                "reportInterval": self.reportInterval,
                "enforcePeriodicBox": True,
            },
            "xtc": {
                "file": "trajectory{}xtc".format(self.runID_str),
                "reportInterval": self.reportInterval,
                "enforcePeriodicBox": True,
            },
            "restart": {
                "reportInterval": 100000,
                "xml": "restart{}xml".format(self.runID_str),
                "chk": "restart{}chk".format(self.runID_str),
            },
        }

        if numberOfSteps is not None:
            self.reporters["screen"]["reportInterval"] = max(1, numberOfSteps // 100)
            self.reporters["screen"]["totalSteps"] = numberOfSteps
            self.reporters["screen"]["progress"] = True
            self.reporters["screen"]["remainingTime"] = True

        self.active_reporters = ["log"]

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

        if simulation is None:
            return

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

            elif rep == "restart":
                if self.reporters[rep]["xml"] is not None:
                    self.addChekpoint(
                        simulation,
                        {
                            "file": self.reporters[rep]["xml"],
                            "reportInterval": self.reporters[rep]["reportInterval"],
                            "writeState": True,
                        },
                    )
                if self.reporters[rep]["chk"] is not None:
                    self.addChekpoint(
                        simulation,
                        {
                            "file": self.reporters[rep]["chk"],
                            "reportInterval": self.reporters[rep]["reportInterval"],
                            "writeState": False,
                        },
                    )
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

    def addChekpoint(self, simulation, config):
        my.pretty_log(
            config,
            title="Checkpoint reporter:",
            indent=1,
            logger="debug",
        )
        simulation.reporters.append(app.CheckpointReporter(**config))
