### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import numpy as np
import openmm as mm
import openmm.unit as unit
from openmm.unit import MOLAR_GAS_CONSTANT_R

import new_openmm_wrapper as my


def computeEnergyFromContext(
    c, outputUnit=unit.kilojoule_per_mole, quiet=False, logger="INFO"
):
    if not quiet:
        my.pretty_log(title="Single point energy calculation", sep=True, logger=logger)

    total_energy = (
        c.getState(getEnergy=True).getPotentialEnergy().value_in_unit(outputUnit)
    )
    if quiet:
        return total_energy

    # Set force groups
    energy_breakdown = {}
    n = 1
    for i, f in enumerate(c.getSystem().getForces()):
        if isinstance(f, mm.CMMotionRemover):
            continue
        f.setForceGroup(i + n)
        energy_breakdown[f.getName()] = {i + n}
        if isinstance(f, mm.NonbondedForce):
            n += 1
            f.setReciprocalSpaceForceGroup(i + n)
            energy_breakdown["ReciprocalSpaceForce"] = {i + n}
            energy_breakdown["Total electrostatics"] = {i + n - 1, i + n}

    my.pretty_log(
        title="Creating force groups:",
        data=energy_breakdown,
        align_width=0,
        logger="debug",
    )
    c.reinitialize(preserveState=True)

    parameters = c.getParameters()
    if len(parameters) > 0:
        my.pretty_log(
            dict(parameters),
            title="Global parameters in context:",
            indent=1,
            logger=logger,
        )

    for name, idx in energy_breakdown.items():
        e = (
            c.getState(getEnergy=True, groups=idx)
            .getPotentialEnergy()
            .value_in_unit(outputUnit)
        )
        energy_breakdown[name] = e
    energy_breakdown["Total Potential Energy"] = total_energy

    symbol = outputUnit.get_symbol()
    log_str = {}
    for k, v in energy_breakdown.items():
        log_str[f"{k} ({symbol})"] = f"{v:15.4f}"

    my.pretty_log("Energy breakdown:", indent=1, logger=logger)
    my.pretty_log(log_str, indent=2, logger=logger)

    return total_energy

    # if eKin:
    #     ndof = my.get_number_of_degrees_of_freedom(s)
    #     k = c.getState(getEnergy=True).getKineticEnergy()
    #     if k.value_in_unit(outputUnit) > 1e-3:
    #         logger.info(
    #             "  {:40s} = {} ".format(
    #                 "Kinetic energy (" + symbol + ")", k.value_in_unit(outputUnit)
    #             )
    #         )
    #         temp = 2 * k / ndof / MOLAR_GAS_CONSTANT_R
    #         logger.info("  {:40s} = {} ".format("Temperature (K)", temp.value_in_unit(unit.kelvin)))
