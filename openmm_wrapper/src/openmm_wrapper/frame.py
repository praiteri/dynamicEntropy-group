# -*- coding: utf-8 -*-
"""This module defines a class for handling trajectory frames."""

__all__ = ["Frame"]


class Frame(object):
    """A class for storing trajectory frame coordinates and provide methods
    acting on them."""

    def __init__(self, index, natoms, coords, unitcell=None, velocs=None):

        self._index = index
        self._natoms = natoms
        self._coords = coords
        self._velocs = velocs
        self._unitcell = unitcell

    def numAtoms(self):
        """Returns number of atoms."""
        return self._natoms

    def getPositions(self):
        return self._coords

    def getCell(self):
        return self._unitcell
