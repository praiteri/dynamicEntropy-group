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
from openmm.unit import MOLAR_GAS_CONSTANT_R

import openmm_wrapper as my


def computeEnergyFromContext(
    c, outputUnit=unit.kilojoule_per_mole, header=None, eKin=False
):
    logger = logging.getLogger("dynamicEntropy")
    if header is None:
        logger.critical("#--- Single point energy ------------------#")
    else:
        logger.info(header)

    logger.info("  Context parameters ...")
    for p in c.getParameters():
        logger.info("    {:38s} = {} ".format(p, c.getParameter(p)))

    s = c.getSystem()
    e = c.getState(getEnergy=True).getPotentialEnergy().value_in_unit(outputUnit)

    symbol = outputUnit.get_symbol()
    logger.info("  {:40s} = {} ".format("Total potential energy (" + symbol + ")", e))

    for force in s.getForces():
        if isinstance(force, mm.CMMotionRemover):
            continue
        e = (
            c.getState(getEnergy=True, groups={force.getForceGroup()})
            .getPotentialEnergy()
            .value_in_unit(outputUnit)
        )
        logger.info("  {:40s} = {} ".format(force.getName() + " (" + symbol + ")", e))

    if eKin:
        ndof = my.getNumberOfDegreesOfFreedom(s)
        k = c.getState(getEnergy=True).getKineticEnergy()
        if k.value_in_unit(outputUnit) > 1e-3:
            logger.info(
                "  {:40s} = {} ".format(
                    "Kinetic energy (" + symbol + ")", k.value_in_unit(outputUnit)
                )
            )
            temp = 2 * k / ndof / MOLAR_GAS_CONSTANT_R
            logger.info(
                "  {:40s} = {} ".format(
                    "Temperature (K)", temp.value_in_unit(unit.kelvin)
                )
            )

    return
