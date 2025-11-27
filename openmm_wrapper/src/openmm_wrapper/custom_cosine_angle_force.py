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

import openmm as mm
import openmm.app as app
import openmm.unit as unit


def check_value(key, value, list_of_dicts):
    matches = []
    for index, dictionary in enumerate(list_of_dicts):
        if key in dictionary and dictionary[key] == value:
            matches.append(index)
    return matches


def intersection_with_duplicates(type1, type2):
    """
    Compute the intersection of two lists with duplicates.
    """
    from collections import Counter

    counter1 = Counter(type1)
    counter2 = Counter(type2)

    # Compute the intersection by taking the minimum count for each element
    intersection_result = {
        element: min(counter1[element], counter2[element])
        for element in set(type1) & set(type2)
    }

    # Convert the result back to a list or tuple if needed
    result_list = [
        element for element, count in intersection_result.items() for _ in range(count)
    ]

    return result_list


def customCosineAngleForce(system, setup, forcefield, topology, customCosineAngle):
    import itertools
    from collections import Counter

    expression = "0.5*k*cos(n*theta-theta0)^2"
    force = mm.CustomAngleForce(expression)

    force.addPerAngleParameter("k")
    force.addPerAngleParameter("n")
    force.addPerAngleParameter("theta0")

    listOfAtoms = setup.getListOfAtoms()
    bondedToAtom = forcefield._buildBondedToAtomList(topology)

    for a in listOfAtoms:
        for index in check_value("pivot", a["type"], customCosineAngle):
            angle = customCosineAngle[index]
            angleTypes = [
                angle[t]
                for t in ["type1", "type2"]
                if t in angle and "*" not in angle[t]
            ]
            listOfBondedPairs = list(
                itertools.combinations(bondedToAtom[a["index"]], 2)
            )
            listOfBondedTypes = [
                [listOfAtoms[i[0]]["type"], listOfAtoms[i[1]]["type"]]
                for i in listOfBondedPairs
            ]
            prm = [
                angle["params"]["k"],
                angle["params"]["n"],
                angle["params"]["theta0"],
            ]
            for t, i in zip(listOfBondedTypes, listOfBondedPairs):
                if Counter(t) == Counter(angleTypes):
                    force.addAngle(i[0], a["index"], i[1], prm)

    force.setName("cosineAngleForce")
    system.addForce(force)
