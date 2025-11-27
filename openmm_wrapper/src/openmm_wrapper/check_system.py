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
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit


def checkSystem(system, topology, forcefield, pos):
    """
    Prints information about the system for identifying simple issues with the input.
    """

    logger = logging.getLogger("dynamicEntropy")
    try:
        nonbonded = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][
            0
        ]

        def getCharge(i):
            return nonbonded.getParticleParameters(i)[0].value_in_unit(
                unit.elementary_charge
            )

    except:

        def getCharge(i):
            return 0

        pass

    try:
        nonbonded = [
            f for f in system.getForces() if isinstance(f, mm.AmoebaMultipoleForce)
        ][0]

        def getCharge(i):
            return nonbonded.getMultipoleParameters(i)[0].value_in_unit(
                unit.elementary_charge
            )

    except:
        pass

    logger.critical("System details ...")

    try:
        cell = topology.getPeriodicBoxVectors().value_in_unit(unit.nanometer)
        logger.info(
            "  {:40s} = {:-7.4f} {:-7.4f} {:-7.4f}".format(
                "Periodic cell vector a (nm)", *cell[0]
            )
        )
        logger.info(
            "  {:40s} = {:-7.4f} {:-7.4f} {:-7.4f}".format(
                "Periodic cell vector b (nm)", *cell[1]
            )
        )
        logger.info(
            "  {:40s} = {:-7.4f} {:-7.4f} {:-7.4f}".format(
                "Periodic cell vector c (nm)", *cell[2]
            )
        )
        logger.info("  {:40s} = {}".format("Number of atoms", topology.getNumAtoms()))
    except:
        pass

    residueNames = []
    listOfResidues = []
    for res, r in zip(forcefield.getMatchingTemplates(topology), topology.residues()):
        if res.name not in residueNames:
            residueNames.append(res.name)
            listOfResidues.append([res, r, 1])
        else:
            idx = residueNames.index(res.name)
            listOfResidues[idx][-1] += 1

    logger.info("  {:40s} = {}".format("Number of residue types", len(listOfResidues)))

    listOfAtoms = []
    totalCharge = 0.0
    for res in listOfResidues:
        logger.info(
            "    {:38s} = {}".format(
                "Number of residues of type " + res[0].name, res[2]
            )
        )
        for ap, af, r in zip(res[1]._atoms, res[0].atoms, pos):
            charge = getCharge(ap.index)
            totalCharge += charge * res[2]
            listOfAtoms.append([ap.name, charge, af.type, af.name])
            logger.debug(
                "Atom {:3s} - charge {:+10.7f} type   = {}/{}".format(*listOfAtoms[-1])
            )

    logger.info("  {:40s} = {:<10.3f}".format("Total charge", totalCharge))

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
