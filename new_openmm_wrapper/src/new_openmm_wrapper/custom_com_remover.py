### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import numpy as np
import logging
import openmm as mm
import openmm.app as app
import new_openmm_wrapper as my

import openmm.unit as unit


class customCOMRemover(object):
    def __init__(self, settings, topology):
        self.logger = logging.getLogger("dynamicEntropy")
        self.logger.debug("Initializing custom COM remover")
        self.logger.debug(f"Settings: {settings}")

        self._reportInterval = int(settings["reportInterval"])

        self._atom_indices = [
            my.get_atoms_from_topology({"species": g}, topology) for g in settings["groups"]
        ]

        self._atoms_masses = np.array(
            [atom.element.mass.value_in_unit(unit.dalton) for atom in topology.atoms()]
        )

        self._group_mass = np.array(
            [self._atoms_masses[indices].sum() for indices in self._atom_indices]
        )

        self.logger.info(
            f"Custom COM remover will act on {len(self._atom_indices)} groups of atoms, every {self._reportInterval} steps"
        )

    def __del__(self):
        return

    def describeNextReport(self, simulation):
        """
        # Returns a five element tuple.
        # The first element is the number of steps until the next report.
        # The remaining elements specify whether that report will require
        # positions, velocities, forces, and energies respectively
        """
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, True, False, False, None)

    def report(self, simulation, state):
        # stime = state.getTime().value_in_unit(unit.picosecond)

        # forceGroup = set()
        state = simulation.context.getState(getVelocities=True)

        # positions = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
        velocities = state.getVelocities(asNumpy=True).value_in_unit(
            unit.nanometer / unit.picosecond
        )

        for i, atom_indices in enumerate(self._atom_indices):
            masses = self._atoms_masses[atom_indices]
            group_momentum = np.dot(masses, velocities[atom_indices])
            com_velocity = group_momentum / self._group_mass[i]
            velocities[atom_indices] -= com_velocity

        velocities = unit.Quantity(velocities, unit.nanometer / unit.picosecond)

        simulation.context.setVelocities(velocities)

        return
