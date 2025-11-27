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

import logging
import openmm as mm
import openmm.unit as unit
import openmm_wrapper as my


def integrator(system, config):
    logger = logging.getLogger("dynamicEntropy")

    if config["thermostat"].upper() in {"LANG", "VLANG"}:
        if config["thermostatParameter"].unit == unit.picosecond:
            config["thermostatParameter"] = config["thermostatParameter"]._value * (
                1 / unit.picosecond
            )

    argsThermo = [
        config["temperature"],
        config["thermostatParameter"],
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

            # elif config["barostat"].upper() in ["MEMBRANE"]:
            #     system.addForce(mm.MonteCarloMembraneBarostat(
            #         config["pressure"], 0.0, config["temperature"],
            #         1 , 0,
            #         config["barostatUpdate"])
            #                     )

            else:
                raise Exception("Unknown barostat")

        if config["thermostat"].upper() == "LANG":
            integrator = mm.LangevinMiddleIntegrator(*argsThermo)

        elif config["thermostat"].upper() == "VLANG":
            integrator = mm.VariableLangevinIntegrator(*argsThermo)
            logger.warning(
                "Running Variable LANG with maxStepSize of %s. (0 -> No limit)",
                integrator.getMaximumStepSize(),
            )

        elif config["thermostat"].upper() == "CSVR":
            integrator = my.CSVRIntegrator(*argsThermo, system, max_d=maxDisplacement)
            if maxDisplacement is not None:
                logger.warning(
                    "Running CSVR with a maxDisplacement of %s", maxDisplacement
                )

        elif config["thermostat"].upper() == "NH":
            integrator = mm.NoseHooverIntegrator(*argsThermo, 3)

        # elif config["md"]["integrator"].upper() == "AND":
        #     integrator = mm.VerletIntegrator(config["timestep"])
        #     system.addForce(mm.AndersenThermostat(*argsThermo[0:2]))

        else:
            raise Exception(
                "Unknown integrator/thermostat ({})".format(config["thermostat"])
            )

    return integrator
