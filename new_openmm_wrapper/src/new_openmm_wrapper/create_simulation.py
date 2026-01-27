### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import openmm as mm
import openmm.app as app
import openmm.unit as unit
import new_openmm_wrapper as my


def createSimulation(topology, system, integrator, platform, properties):
    simulation = app.Simulation(topology, system, integrator, platform, properties)

    try:
        nonbonded = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][0]
    except Exception:
        nonbonded = [f for f in system.getForces() if isinstance(f, mm.AmoebaMultipoleForce)][0]

    nbm = nonbonded.getNonbondedMethod()
    nbm_list = {
        0: "NoCutoff",
        1: "CutoffNonPeriodic",
        2: "CutoffPeriodic",
        3: "Ewald",
        4: "PME",
        5: "LJPME",
    }

    my.pretty_log(title="Creating simulation object", sep=True)
    my.pretty_log({"Nonbonded method": nbm_list[nbm]}, indent=1)
    my.pretty_log({"Global cutoff distance": nonbonded.getCutoffDistance()}, indent=1)
    try:
        my.pretty_log({"Use switching distance": nonbonded.getUseSwitchingFunction()}, indent=1)
        my.pretty_log({"Switching distance": nonbonded.getSwitchingDistance()}, indent=1)
    except Exception:
        pass

    if nbm in [4, 5]:
        my.pretty_log(
            {"PME parameters": nonbonded.getPMEParametersInContext(simulation.context)}, indent=1
        )

    if nbm in [3, 4, 5]:
        my.pretty_log({"Ewald tolerance": nonbonded.getEwaldErrorTolerance()}, indent=1)

    if nbm in [2]:
        my.pretty_log(
            {"Ewald tolerance": nonbonded.getReactionFieldDielectric()}, logger="warning", indent=1
        )

    return simulation
