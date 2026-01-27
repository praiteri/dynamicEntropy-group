### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import openmm as mm
import openmm.app as app
import openmm.unit as unit
import new_openmm_wrapper as my


def checkSystem(system, topology, forcefield, pos):
    """
    Prints information about the system for identifying simple issues with the input.
    """

    def getCharge(i):
        """Get particle charge from either NonbondedForce or AmoebaMultipoleForce."""
        try:
            nonbonded = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][0]
            return nonbonded.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
        except Exception:
            pass

        try:
            multipole = [f for f in system.getForces() if isinstance(f, mm.AmoebaMultipoleForce)][0]
            return multipole.getMultipoleParameters(i)[0].value_in_unit(unit.elementary_charge)
        except Exception:
            return 0

    my.pretty_log(title="System details", sep=True)
    try:
        cell = topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer)
        a_str = f"{cell[0][0]:-7.4f} {cell[0][1]:-7.4f} {cell[0][2]:-7.4f}"
        b_str = f"{cell[1][0]:-7.4f} {cell[1][1]:-7.4f} {cell[1][2]:-7.4f}"
        c_str = f"{cell[2][0]:-7.4f} {cell[2][1]:-7.4f} {cell[2][2]:-7.4f}"
        my.pretty_log({"Periodic cell vector a (nm)": a_str}, indent=1)
        my.pretty_log({"Periodic cell vector b (nm)": b_str}, indent=1)
        my.pretty_log({"Periodic cell vector c (nm)": c_str}, indent=1)
    except Exception:
        pass

    my.pretty_log({"Number of atoms": topology.getNumAtoms()}, indent=1)

    residueNames = []
    listOfResidues = []
    for res, r in zip(forcefield.getMatchingTemplates(topology), topology.residues()):
        if res.name not in residueNames:
            residueNames.append(res.name)
            listOfResidues.append([res, r, 1])
        else:
            idx = residueNames.index(res.name)
            listOfResidues[idx][-1] += 1

    my.pretty_log({"Number of residue types": len(listOfResidues)}, indent=1)

    listOfAtoms = []
    totalCharge = 0.0
    for res in listOfResidues:
        my.pretty_log({f"Number of residues of type {res[0].name}": res[2]}, indent=2)
        for ap, af, r in zip(res[1]._atoms, res[0].atoms, pos):
            charge = getCharge(ap.index)
            totalCharge += charge * res[2]
            listOfAtoms.append([ap.name, charge, af.type, af.name])
            log_str = f"Atom {ap.name:4s} : chanrge = {charge:+10.7f}, type = {af.type}/{af.name}"
            my.pretty_log(log_str, indent=3, logger="debug")

    my.pretty_log({"Total charge": round(totalCharge, 3)}, indent=1)

    # q = 0.0
    # x = mm.Vec3(0,0,0) * unit.nanometers / 2
    # dipole = mm.Vec3(0,0,0) * unit.nanometers
    # # for i in range(system.getNumParticles()):
    # for i in range(3):
    #     q += getCharge(i)
    #     dipole += getCharge(i) * (pos[i] - x)
    # q *= unit.elementary_charge
    # logger.info("  {:40s} = {:<10.3f}".format("Total charge (e)",q.value_in_unit(unit.elementary_charge)))
    # dipole *= unit.elementary_charge
    # logger.info("  {:40s} = {:<10.3f} {:<10.3f} {:<10.3f}".format("Total dipole (Debye)",*dipole.value_in_unit(unit.debye)))
