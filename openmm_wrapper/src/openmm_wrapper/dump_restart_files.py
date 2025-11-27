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
import os


def dumpRestartFiles(simulation, fileName, XML=True, PDB=True):
    """
    Write restart file in XML format and final coordinates in PDB
    """
    if XML:
        simulation.saveState(fileName)

    if PDB:
        state = simulation.context.getState(getPositions=True)
        pos = state.getPositions()
        base = os.path.splitext(fileName)
        simulation.topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())

        ff = base[0] + ".pdb"
        app.PDBFile.writeFile(simulation.topology, pos, open(ff, "w"))
