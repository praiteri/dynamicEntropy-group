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


def customDistanceImproperForce(sys, setup, ffx, topology, customDistanceImproper):
    """
    Improper torsion force defined as the distance of the central atom
    from the plane formed by the other three (U = Ad^2 + Bd^4).
    Written by Blake I. Armstrong - Curtin Univeristy
    """
    expression = "k2*dst^2+k4*dst^4"
    expression += "; dst=nui*ci+nuj*cj+nuk*ck"
    expression += "; nui=(ni)/(sqrt(ni^2+nj^2+nk^2))"
    expression += "; nuj=(nj)/(sqrt(ni^2+nj^2+nk^2))"
    expression += "; nuk=(nk)/(sqrt(ni^2+nj^2+nk^2))"
    expression += "; ni=((aj*bk)-(ak*bj))"
    expression += "; nj=-((ai*bk)-(ak*bi))"
    expression += "; nk=((ai*bj)-(aj*bi))"
    expression += "; ai=x4-x2"
    expression += "; aj=y4-y2"
    expression += "; ak=z4-z2"
    expression += "; bi=x4-x3"
    expression += "; bj=y4-y3"
    expression += "; bk=z4-z3"
    expression += "; ci=x4-x1"
    expression += "; cj=y4-y1"
    expression += "; ck=z4-z1"

    tor = mm.CustomCompoundBondForce(4, expression)
    tor.addPerBondParameter("k2")
    tor.addPerBondParameter("k4")

    listOfAtoms = setup.getListOfAtoms()
    bondedToAtom = ffx._buildBondedToAtomList(topology)

    addForce = False
    # With the loops in this order only the first match is used
    for a in listOfAtoms:
        if len(bondedToAtom[a["index"]]) != 3:
            continue

        for index in check_value("centre", a["type"], customDistanceImproper):
            lst = [a["index"]] + bondedToAtom[a["index"]]
            impr = customDistanceImproper[index]
            bondedTypes = [listOfAtoms[i]["type"] for i in lst[1:]]
            improperTypes = [
                impr[t]
                for t in ["type1", "type2", "type3"]
                if t in impr and "*" not in impr[t]
            ]
            match = intersection_with_duplicates(bondedTypes, improperTypes)

            if len(match) != len(improperTypes):
                continue

            addForce = True
            k2 = impr["params"]["k2"]
            k4 = impr["params"]["k4"]
            if "conv" in impr:
                raise Exception(
                    "conv has been removed, add conversion to the parameters"
                )

            tor.addBond(lst, [k2, k4])

    if addForce:
        tor.setName("DistanceImproper")
        sys.addForce(tor)
    else:
        del tor
