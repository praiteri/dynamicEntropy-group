### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import logging
import openmm as mm
import openmm.unit as unit
import new_openmm_wrapper as my


def createIntegrator(system, config):
    argsThermo = [
        config["temperature"],
        config["thermostatParameter"]._value,
        config["timestep"],
    ]
    maxDisplacement = config.get("maxDisplacement", None)

    if config["ensemble"].upper() == "NVE":
        integrator = mm.VerletIntegrator(config["timestep"])

    else:
        if config["ensemble"].upper() == "NPT":
            if config["barostat"].upper() == "TRI":
                # logger.critical("#--- Using triclinic barostat -------#")
                system.addForce(
                    mm.MonteCarloFlexibleBarostat(
                        config["pressure"],
                        config["temperature"],
                        config["barostatUpdate"],
                    )
                )

            elif config["barostat"].upper() == "ISO":
                # logger.critical("#--- Using isotropic barostat -------------#")
                system.addForce(
                    mm.MonteCarloBarostat(
                        config["pressure"],
                        config["temperature"],
                        config["barostatUpdate"],
                    )
                )

            elif config["barostat"].upper() in [
                "ORTHO",
                "X",
                "Y",
                "Z",
                "XY",
                "XZ",
                "YZ",
            ]:
                # logger.critical("#--- Using orthorhombic barostat ----------#")
                press = (config["pressure"],) * 3
                if config["barostat"].upper() == "ORTHO":
                    flx = [True, True, True]
                elif config["barostat"].upper() == "X":
                    flx = [True, False, False]
                elif config["barostat"].upper() == "Y":
                    flx = [False, True, False]
                elif config["barostat"].upper() == "Z":
                    flx = [False, False, True]
                elif config["barostat"].upper() == "XY":
                    flx = [True, True, False]
                elif config["barostat"].upper() == "XZ":
                    flx = [True, False, True]
                elif config["barostat"].upper() == "YZ":
                    flx = [False, True, True]
                system.addForce(
                    mm.MonteCarloAnisotropicBarostat(
                        press, config["temperature"], *flx, config["barostatUpdate"]
                    )
                )

            else:
                raise Exception("Unknown barostat")

        if config["thermostat"].upper() == "LANG":
            integrator = mm.LangevinMiddleIntegrator(*argsThermo)

        elif config["thermostat"].upper() == "VLANG":
            integrator = mm.VariableLangevinIntegrator(*argsThermo)
            my.pretty_log(
                {"Running Variable LANG with maxStepSize": integrator.getMaximumStepSize()},
                logger="warning",
            )

        elif config["thermostat"].upper() == "CSVR":
            integrator = my.CSVRIntegrator(*argsThermo, system, max_d=maxDisplacement)
            if maxDisplacement is not None:
                my.pretty_log(
                    {"Running CSVR with maxDisplacement": maxDisplacement},
                    logger="warning",
                )

        elif config["thermostat"].upper() == "NH":
            integrator = mm.NoseHooverIntegrator(*argsThermo, 3)

        else:
            raise Exception("Unknown integrator/thermostat ({})".format(config["thermostat"]))

    log_str = {"Timestep": config["timestep"]}
    if config["ensemble"].upper() in ["NVT", "NPT"]:
        log_str["Temperature"] = config["temperature"]
        log_str["Thermostat"] = config["thermostat"].upper()
        log_str["Thermostat Parameter"] = config["thermostatParameter"]
    if config["ensemble"].upper() in ["NPT"]:
        log_str["Pressure"] = config["pressure"]
        log_str["Barostst"] = config["barostat"].upper()
        log_str["Barostat Update Freq"] = config["barostatUpdate"]
    my.pretty_log(log_str, title="Creating integrator")

    if config["thermostat"].upper() in {"LANG", "VLANG"}:
        if config["thermostatParameter"].unit == unit.picosecond:
            config["thermostatParameter"] = config["thermostatParameter"]._value * (
                1 / unit.picosecond
            )

    _ = my.get_number_of_degrees_of_freedom(system, "debug")

    return integrator
