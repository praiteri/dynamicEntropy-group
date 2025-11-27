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

import openmm_wrapper as my

import os
import sys
import numpy as np
import logging
from pprint import *


class CatalyticReactionReporter(object):
    def __init__(self, config):
        """
        .
        """
        logger = logging.getLogger("dynamicEntropy")
        logger.info("osmotic pressure calculation ...")

        self._surface_mask = None

        self._bonded_pairs = set()
        self._bonded_atoms = set()

        self._localCounter = 0

        self._computeInterval = 1000
        if "computeInterval" in config:
            self._computeInterval = config["computeInterval"]

        self.substrate = config.get("substrate", None)
        if self.substrate is None:
            raise ValueError("CatalyticReactionReporter: substrate not specified")

        self._reactants = config.get("reactants", None)
        if self._reactants is None:
            raise ValueError("CatalyticReactionReporter: reactants not specified")
        if not isinstance(self._reactants, list):
            self._reactants = self._reactants.split(",")
        self._reactants = [r.strip() for r in self._reactants]

        logger.info("CATALYTIC REACTION: --------------------------#")
        logger.debug(pformat(config))
        logger.debug(f"Compute interval: {self._computeInterval} steps")
        logger.info(f"Substrate: {self.substrate}")
        logger.info(f"Reactants: {self._reactants}")
        logger.info("CATALYTIC REACTION: --------------------------#")

    def __del__(self):
        try:
            self._out.close()
        except:
            pass

    def describeNextReport(self, simulation):
        """
        # Returns a five element tuple.
        # The first element is the number of steps until the next report.
        # The remaining elements specify whether that report will require
        # positions, velocities, forces, and energies respectively
        """
        steps = self._computeInterval - simulation.currentStep % self._computeInterval
        return (steps, True, False, False, False)

    def report(self, simulation, state):
        self._localCounter += 1

        state = simulation.context.getState(getPositions=True)
        positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)

        # Create boolean mask for non-Au atoms
        if self._surface_mask is None:
            atom_names = np.array([a.name for a in simulation.topology.atoms()])
            self._surface_mask = np.isin(atom_names, self.substrate)

        # Apply mask and z-coordinate filter in one go
        # Create combined condition
        z_layer = 2.0  # nm
        combined_mask = self._surface_mask & (positions[:, 2] < z_layer)

        # Get indices and positions
        final_indices = np.where(combined_mask)[0]
        r_surf = positions[combined_mask]

        box_lengths = np.array(
            state.getPeriodicBoxVectors().value_in_unit(unit.nanometer)
        ).diagonal()

        neighbour_distances, (i_idx, j_idx) = self.pairwise_distances_pbc_upper_triangle(
            r_surf, box_lengths, final_indices
        )

        nrand = np.random.uniform(0, 1, len(neighbour_distances))

        bonded_pairs_set = self._bonded_pairs | {(j, i) for i, j in self._bonded_pairs}
        bonded_atoms_set = set(self._bonded_atoms)

        bonds = [
            (int(i), int(j), d, n)
            for i, j, d, n in zip(i_idx, j_idx, neighbour_distances, nrand)
            if d < 0.25
            and (int(i), int(j)) not in bonded_pairs_set
            and int(i) not in bonded_atoms_set
            and int(j) not in bonded_atoms_set
        ]

        print(f"Step {simulation.currentStep}: Found {len(bonds)} new bonds")

        system = simulation.context.getSystem()

        for force in system.getForces():
            if force.__class__.__name__ == "HarmonicBondForce":
                bondForce = force
            # if force.__class__.__name__ == "CustomNonbondedForce":
            #     nonBondedForce = force

        # bonds = []
        # for b in simulation.topology.bonds():
        #     self.listOfBonds.append((b.atom1.index, b.atom2.index))
        # exclusions = app.forcefield._findExclusions(bonds, 3, system.getNumParticles())

        if len(bonds) == 0:
            return

        # n = bondForce.getNumBonds()
        # if n > 0:
        #     print(f"Number of bonds: {n}")
        #     print("Existing bonds:")
        #     for b in self._bonded_pairs:
        #         print(b)
        # print(f"Bonded atom: {self._bonded_atoms}")

        for i, j, d, n in bonds:
            # if n > 0.5:
            #     continue

            # Bonds cannot be added/removed during a simulation
            # We can only modify the parameters of existing bonds with setBondParameters()
            # and then call updateParametersInContext()
            # To add bonds we can use addBond() but then we have to recreate the Context
            # context.reinitialize(self, preserveState=True), which should preserve positions and velocities

            self._bonded_pairs.add((i, j))
            self._bonded_atoms.add(i)
            self._bonded_atoms.add(j)
            bondForce.addBond(int(i), int(j), 0.2, 40000.0)

            for force in system.getForces():
                try:
                    force.addExclusion(int(i), int(j))
                except Exception:
                    pass

                try:
                    force.addException(int(i), int(j), 0, 0, 1)
                except Exception:
                    pass

            print(f"Added bond between {i} and {j}, distance: {d:.3f} nm")
            print(f"Positions: {positions[i]}, {positions[j]}")
            simulation.context.reinitialize(preserveState=True)
            # quit()

    def pairwise_distances_pbc_upper_triangle(self, r_surf, box_lengths, map):
        """
        Compute only upper triangle of distance matrix (since it's symmetric).
        Returns flattened array of unique pairwise distances.

        Parameters:
        -----------
        r_surf : numpy.ndarray, shape (N, 3)
            Atomic coordinates where N is the number of atoms
        box_lengths : numpy.ndarray, shape (3,)
            Box dimensions [Lx, Ly, Lz] for the orthogonal cell

        Returns:
        --------
        distances : numpy.ndarray, shape (N*(N-1)/2,)
            Flattened array of unique pairwise distances
        indices : tuple of arrays
            Row and column indices corresponding to the distances
        """
        N = r_surf.shape[0]

        # Get upper triangle indices (excluding diagonal)
        i_idx, j_idx = np.triu_indices(N, k=1)

        # Get coordinates for pairs
        r_i = r_surf[i_idx]  # Shape (n_pairs, 3)
        r_j = r_surf[j_idx]  # Shape (n_pairs, 3)

        # Compute differences
        dr = r_i - r_j

        # Apply PBC
        dr = dr - box_lengths * np.round(dr / box_lengths)

        # Compute distances
        distances = np.sqrt(np.sum(dr**2, axis=1))

        return distances, (map[i_idx], map[j_idx])
