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

import openmm.unit as unit


def addElectronVoltToUnit():
    r0 = (1 * unit.elementary_charge * unit.volt).value_in_unit(unit.joule)
    r1 = (
        1 * unit.elementary_charge * unit.volt * unit.AVOGADRO_CONSTANT_NA
    ).value_in_unit(unit.kilojoule_per_mole)

    ev_base_unit = unit.ScaledUnit(r0, unit.kilojoule_per_mole, "electron_volt", "eV")
    unit.ev_base_unit = ev_base_unit
    unit.ev = unit.electron_volt = unit.Unit({ev_base_unit: 1.0})

    md_ev_base_unit = unit.ScaledUnit(
        r1, unit.kilojoule_per_mole, "electron_volt", "eV"
    )
    unit.md_ev_base_unit = md_ev_base_unit
    unit.md_ev = unit.md_electron_volt = unit.Unit({md_ev_base_unit: 1.0})
