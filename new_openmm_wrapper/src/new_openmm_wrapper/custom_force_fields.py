### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import logging


def createCustomForceFields(customFunctionsFromFile=None):
    logger = logging.getLogger("dynamicEntropy")

    customForceFields = {
        "lj": {
            "type": "lj",
            "params": ["e", "s"],
            "init": {"e": 0, "s": 1},
            "expression": ("4*e*((s/r)^12-(s/r)^6);e=e(type1,type2); s=s(type1, type2)"),
            "fep": (
                "{var}^2*4*e*(1./(alpha+(r/s)^6)^2-1/(alpha+(r/s)^6));"
                "alpha=(1-{var})^2;"
                "e=e(type1,type2);"
                "s=s(type1, type2)"
            ),
            "fep2": ("{var}*4*e*((s/r)^12-(s/r)^6);e=e(type1,type2); s=s(type1, type2)"),
        },
        "buck": {
            "type": "buck",
            "params": ["A", "rho", "C"],
            "init": {"A": 0, "rho": 1, "C": 0},
            "expression": (
                "A*exp(-r/rho)-C/r^6; A=A(type1, type2); rho=rho(type1, type2); C=C(type1, type2)"
            ),
            "fep": (
                "{var}^4*(A*exp(-r/rho) - C/(alpha+r^2)^3);"
                "alpha=(1-{var})^2;"
                "A=A(type1, type2);"
                "rho=rho(type1, type2);"
                "C=C(type1, type2)"
            ),
            "fep2": (
                "{var}*A*exp(-r/rho)-C/r^6; A=A(type1, type2); rho=rho(type1, type2); C=C(type1, type2)"
            ),
        },
        "AB": {
            "type": "AB",
            "params": ["A", "B"],
            "init": {"A": 0, "B": 0},
            "expression": ("A/r^12 - B/r^6; A=A(type1,type2); B=B(type1,type2)"),
            "fep": (
                "{var}^2*(A/(alpha+r^2)^6-B/(alpha+r^2)^3);"
                "alpha=(1-{var})^2;"
                "A=A(type1,type2);"
                "B=B(type1,type2)"
            ),
            "fep2": ("{var}*A/r^12 - B/r^6; A=A(type1,type2); B=B(type1,type2)"),
        },
        "harmonic": {
            "type": "harmonic",
            "params": ["k", "r0"],
            "init": {"k": 0, "r0": 1},
            "expression": ("step(r0-r)*0.5*k*(r-r0)^2); k=k(type1,type2); r0=r0(type1,type2)"),
            "fep": ("{var}*step(r0-r)*0.5*k*(r-r0)^2; k=k(type1,type2); r0=r0(type1,type2)"),
            "fep2": ("{var}*step(r0-r)*0.5*k*(r-r0)^2; k=k(type1,type2); r0=r0(type1,type2)"),
        },
    }

    if customFunctionsFromFile is None:
        return customForceFields

    if not isinstance(customFunctionsFromFile, list):
        raise Exception(
            "customForceFields defined in the XML script must be a list of dictionaries"
        )

    for n in customFunctionsFromFile:
        if n["type"] in customForceFields:
            customForceFields[n["type"]] = n
            logger.info(
                "Custom forcefield type {} replaced with definition in the XML file".format(
                    n["type"]
                )
            )
        else:
            customForceFields[n["type"]] = n
            logger.info("  {:40s} = {}".format("New custom forcefield type added", n["type"]))

    # customForceFields["FEP_coul"] = {"type":"FEP_coul",
    #                         "params":["qq"],
    #                         "init":[0],
    #                         "expression":"(1-varCoul^2)*qq/r; \
    #                                         qq=qq(type1, type2)"
    #                         }

    return customForceFields
