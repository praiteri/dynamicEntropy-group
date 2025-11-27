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
import pprint as pp
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import numpy as np


def OPLS_MixingRules(setup, system):
    import numpy as np
    import openmm as mm

    # Add force with the correct combination rules
    forces = {
        system.getForce(index).getName(): system.getForce(index)
        for index in range(system.getNumForces())
    }

    if "OPLS" in forces:
        return

    scaleFactor = [(0.00, 0.00), (0.00, 0.00), (0.50, 0.50), (1.00, 1.00)]

    if "OPLS_Scale14" in setup.config["forcefield"]:
        scaleFactor[2] = setup.config["forcefield"]["OPLS_Scale14"]
    if "OPLS_Scale15" in setup.config["forcefield"]:
        scaleFactor[3] = setup.config["forcefield"]["OPLS_Scale15"]

    if any(num < 1 for num in scaleFactor[3]):
        setup.setExclusionsList(system, nbonds=4)

    nonbonded_force = forces["NonbondedForce"]

    lorentz = mm.CustomNonbondedForce(
        "4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)"
    )
    lorentz.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
    lorentz.addPerParticleParameter("sigma")
    lorentz.addPerParticleParameter("epsilon")
    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
    lorentz.setName("OPLS")
    system.addForce(lorentz)

    # Set the correct 1-4 interaction parameters
    # These are stored in the original nonbonded force
    # and are excluded from the new force
    exclusionList = setup.getExclusionsList()
    for l in exclusionList:
        lorentz.addExclusion(*l[0:2])
    # Fix the 1-2, 1-3 and 1-4 interactions

    for i, ex in enumerate([x for x in exclusionList if x[2] <= 3]):
        p1, p2, idx = ex
        q1, s1, e1 = nonbonded_force.getParticleParameters(p1)
        q2, s2, e2 = nonbonded_force.getParticleParameters(p2)
        nonbonded_force.setExceptionParameters(
            i,
            p1,
            p2,
            q1 * q2 * scaleFactor[idx - 1][0],
            np.sqrt(s1 * s2),
            np.sqrt(e1 * e2) * scaleFactor[idx - 1][1],
        )

    # Add any 1-5 interactions
    for i, ex in enumerate([x for x in exclusionList if x[2] == 4]):
        p1, p2, idx = ex
        q1, s1, e1 = nonbonded_force.getParticleParameters(p1)
        q2, s2, e2 = nonbonded_force.getParticleParameters(p2)
        nonbonded_force.addException(
            p1,
            p2,
            q1 * q2 * scaleFactor[idx - 1][0],
            np.sqrt(s1 * s2),
            np.sqrt(e1 * e2) * scaleFactor[idx - 1][1],
        )

    # Set the epsilon of the existing NonBondedForce to zero
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        lorentz.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(index, charge, sigma, epsilon * 0)

    return
