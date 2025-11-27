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
import openmm.app as app
import openmm.unit as unit
import numpy as np
import pprint as pp
import re


def customBondForce(system, setup, interactionsList):

    customBondForces = [
        {
            "type": "morse",
            "params": ["D", "alpha", "r0"],
            "init": {"D": 0, "alpha": 1, "r0": 0},
            "expression": "D*(1-exp(-alpha*(r-r0)))^2",
            "name": "MorseBondForce",
        },
    ]

    # list of atom types
    listOfAtoms = setup.getListOfAtoms()

    # List of bonds in the system
    bonds = setup.getBondsList()

    # List of 1-2 interactions in the system
    list12 = app.forcefield._findExclusions(bonds, 1, system.getNumParticles())

    for cft in customBondForces:
        force = mm.CustomBondForce(cft["expression"])
        for x in cft["params"]:
            _ = force.addPerBondParameter(x)
        force.setName(cft["name"])

        numInt = 0
        for x in list12:
            i1, i2 = x[0:2]
            a1 = listOfAtoms[i1]["type"]
            a2 = listOfAtoms[i2]["type"]

            for iList in interactionsList:
                if (iList["type1"] == a1 and iList["type2"] == a2) or (
                    iList["type1"] == a2 and iList["type2"] == a1
                ):
                    p = [iList["params"][x] for x in cft["params"]]
                    force.addBond(i1, i2, p)
                    numInt += 1

        if numInt > 0:
            system.addForce(force)
        else:
            del force

    return
