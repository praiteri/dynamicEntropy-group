### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import openmm.app as app
import new_openmm_wrapper as my
import os


def getForcefieldFiles(string):
    lst = string.split(",")
    for i, f in enumerate(lst):
        if os.path.exists(f):
            lst[i] = f
        else:
            fr = os.path.join(os.path.dirname(my.__file__), "../forcefields/", f)
            if os.path.exists(fr):
                lst[i] = fr
            else:
                raise Exception("Missing force file file ({})".format(fr))
    return lst


def get_custom_ff_from_file(config):
    customForceFieldsX = None
    # This block should probably move outside this function
    forcefield = my.readForcefield(config)
    if not config["forcefield"]["skipScript"]:
        my.pretty_log("Creating custom forcefields from the XML file")
        for s in forcefield._scripts:
            exec(s, globals())

        for v in ["nonbondedCutoff", "switchDistance"]:
            if v not in globals():
                globals()[v] = config["forcefield"][v]

    if "customForceFieldParameters" in globals():
        # Get customForceFields from globals if it exists, otherwise None
        custom_force_fields = globals().get("customForceFields")

        if custom_force_fields is not None:
            customForceFieldsX = my.createCustomForceFields(
                customFunctionsFromFile=custom_force_fields
            )
        else:
            customForceFieldsX = my.createCustomForceFields()

    my.pretty_log(
        customForceFieldsX,
        title=f"Custom FF extracted from {config['input']['forcefield']}",
        align_width=0,
        logger="debug",
    )

    return customForceFieldsX, globals().get("customForceFieldParameters")


def readForcefield(config):
    fff = getForcefieldFiles(config["input"]["forcefield"])
    forcefield = app.ForceField(*fff)

    for f in forcefield._forces:
        if isinstance(f, app.forcefield.AmoebaMultipoleGenerator):
            config["forcefield"]["isAMOEBA"] = True

    for s in forcefield._scripts:
        for line in s.split("\n"):
            if line.startswith("#@#set"):
                ll = line.split()
                my.locate_and_modify_key(config[ll[1]], ll[2], ll[3])

    return forcefield
