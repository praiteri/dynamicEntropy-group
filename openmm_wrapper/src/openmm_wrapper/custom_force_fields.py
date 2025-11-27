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


def customForceFields(customFunctionsFromFile=None):
    logger = logging.getLogger("dynamicEntropy")

    customForceFields = [
        {
            "type": "lj",
            "params": ["e", "s"],
            "init": {"e": 0, "s": 1},
            "expression": "4*e*( (s/r)^12 - (s/r)^6 ); \
                      e=e(type1, type2); \
                      s=s(type1, type2)",
            "fep": "lambdaVdw^2*4*e*(1./(alpha+(r/s)^6)^2-1/(alpha+(r/s)^6)) ; alpha=(1-lambdaVdw)^2",
            "newexpr": {
                "params": ["e1", "s1"],
                "init": {"e1": 0, "s1": 1},
                "expression": "4*e1*( (s1/r)^12 - (s1/r)^6 ); \
                          e1=e1(type1, type2); \
                          s1=s1(type1, type2)",
            },
        },
        {
            "type": "buck",
            "params": ["A", "rho", "C"],
            "init": {"A": 0, "rho": 1, "C": 0},
            "expression": "A*exp(-r/rho) - C/r^6; \
                      A=A(type1, type2); \
                      rho=rho(type1, type2); \
                      C=C(type1, type2)",
            "fep": "lambdaVdw^4*(A*exp(-r/rho) - C/(alpha+r^2)^3) ; alpha=(1-lambdaVdw)^2",
            "newexpr": {
                "params": ["A1", "rho1", "C1"],
                "init": {"A1": 0, "rho1": 1, "C1": 0},
                "expression": "A1*exp(-r/rho1) - C1/r^6; \
                          A1=A1(type1, type2); \
                          rho1=rho1(type1, type2); \
                          C1=C1(type1, type2)",
            },
        },
        {
            "type": "AB",
            "params": ["A", "B"],
            "init": {"A": 0, "B": 0},
            "expression": "A/r^12 - B/r^6; \
                      A=A(type1, type2); \
                      B=B(type1, type2)",
            "fep": "lambdaVdw^2*(A/(alpha+r^2)^6 - B/(alpha+r^2)^3) ; alpha=(1-lambdaVdw)^2",
            "newexpr": {
                "params": ["A1", "B1"],
                "init": {"A1": 0, "B1": 0},
                "expression": "A1/r^12 - B1/r^6; \
                          A1=A1(type1, type2); \
                          B1=B1(type1, type2)",
            },
        },
        {
            "type": "harmonic",
            "params": ["k", "r0"],
            "init": {"k": 0, "r0": 1},
            "expression": "step(r0-r)*0.5*k*(r-r0)^2; \
                      k=k(type1, type2); \
                      r0=r0(type1, type2)",
            "fep": "lambdaVdw^2*(step(r0-r)*0.5*k*(r-r0)^2)",
            "newexpr": {
                "params": ["k1", "r1"],
                "init": {"k1": 0, "r1": 1},
                "expression": "step(r1-r)*0.5*k1*(r-r1)^2; \
                          k1=k1(type1, type2); \
                          r1=r1(type1, type2)",
            },
        },
    ]

    if customFunctionsFromFile is None:
        return customForceFields

    if not isinstance(customFunctionsFromFile, list):
        raise Exception(
            "customForceFields defined in the XML script must be a list of dictionaries"
        )

    ctypes = [x["type"] for x in customForceFields]
    for n in customFunctionsFromFile:
        if n["type"] in ctypes:
            idx = ctypes.index(n["type"])
            customForceFields[idx] = n
            logger.info(
                "Custom forcefield type {} replaced with definition in the XML file".format(
                    n["type"]
                )
            )
        else:
            customForceFields.append(n)
            logger.info(
                "  {:40s} = {}".format("New custom forcefield type added", n["type"])
            )

    # customForceFields.append({"type":"FEP_coul",
    #                         "params":["qq"],
    #                         "init":[0],
    #                         "expression":"(1-lambdaCoul^2)*qq/r; \
    #                                         qq=qq(type1, type2)"
    #                         })

    return customForceFields
