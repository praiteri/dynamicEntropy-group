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

import openmm.app as app
import openmm_wrapper as my
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


def readForcefield(config):
    fff = getForcefieldFiles(config["input"]["forcefield"])
    forcefield = app.ForceField(*fff)

    for f in forcefield._forces:
        if isinstance(f, app.forcefield.AmoebaMultipoleGenerator):
            config["forcefield"]["isAMOEBA"] = True

    for s in forcefield._scripts:
        for line in s.split("\n"):
            if line.startswith("#@#set"):
                l = line.split()
                my.locate_and_modify_key(config[l[1]], l[2], l[3])

    return forcefield
