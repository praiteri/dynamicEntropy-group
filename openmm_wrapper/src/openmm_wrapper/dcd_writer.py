"""This module defines classes for handling trajectory files in `DCD format`_.

.. _DCD format: http://www.ks.uiuc.edu/Research/namd/2.6/ug/node13.html"""

from struct import calcsize, unpack, pack
from os.path import getsize
import datetime

import numpy as np
from numpy import float32, fromstring

# from dynamicEntropy2.utils.frame import Frame
import logging

PISQUARE = np.pi**2
RECSCALE32BIT = 1
RECSCALE64BIT = 2

now = datetime.datetime.now


class dcdWriter(object):
    def __init__(self, filename, mode="rb", **kwargs):
        if not isinstance(filename, str):
            raise TypeError("filename argument must be a string")
        if not isinstance(mode, str):
            TypeError("mode argument must be string")
        if not mode in ("r", "w", "a", "r+"):
            ValueError("mode string must begin with one of 'r', 'w', 'r+', or " "'a'")

        self._file = None
        self._filename = filename
        if mode in ("a", "r+"):
            self._file = open(filename, "r+b")
            self._file.seek(0)
            mode = "a"
        else:
            if not mode.endswith("b"):
                mode += "b"
            self._file = open(filename, mode)

        self._astype = kwargs.get("astype", None)

        self._nfi = 0
        self._n_csets = 0
        self._n_atoms = 0
        self._unitcell = True

    def __del__(self):
        self._file.close()

    def writeCoord(self, coords, unitcell=None, **kwargs):
        """Write *coords* to a file open in 'a' or 'w' mode.  *coords* may be
        a NUmpy array or a ProDy object that stores or points to coordinate
        data.  Note that all coordinate sets of ProDy object will be written.
        Number of atoms will be determined from the file or based on the size
        of the first coordinate set written.  If *unitcell* is provided for
        the first coordinate set, it will be expected for the following
        coordinate sets as well.  If *coords* is an :class:`~.Atomic` or
        :class:`~.Ensemble` all coordinate sets will be written.

        Following keywords are used when writing the first coordinate set:

        :arg timestep: timestep used for integration, default is 1
        :arg firsttimestep: number of the first timestep, default is 0
        :arg framefreq: number of timesteps between frames, default is 1"""

        if coords.dtype != float32:
            coords = coords.astype(float32)

        n_atoms = coords.shape[-2]
        if self._n_atoms == 0:
            self._n_atoms = n_atoms
        elif self._n_atoms != n_atoms:
            raise ValueError("coords does not have correct number of atoms")

        dcd = self._file
        pack_i_4N = pack("i", self._n_atoms * 4)
        if self._n_csets == 0:
            if unitcell is None:
                self._unitcell = False
            else:
                self._unitcell = True
            timestep = float(kwargs.get("timestep", 1.0))
            first_ts = int(kwargs.get("firsttimestep", 0))
            framefreq = int(kwargs.get("framefreq", 1))
            n_fixed = 0

            pack_i_0 = pack(b"i", 0)
            pack_ix4_0x4 = pack(b"i" * 4, 0, 0, 0, 0)
            pack_i_1 = pack(b"i", 1)
            pack_i_2 = pack(b"i", 2)
            pack_i_4 = pack(b"i", 4)
            pack_i_84 = pack(b"i", 84)
            pack_i_164 = pack(b"i", 164)

            dcd.write(pack_i_84)
            dcd.write(b"CORD")
            dcd.write(pack_i_0)  # 0 Number of frames in file, none written yet
            dcd.write(pack(b"i", first_ts))  # 1 Starting timestep
            dcd.write(pack(b"i", framefreq))  # 2 Timesteps between frames
            dcd.write(pack_i_0)  # 3 Number of timesteps in simulation
            dcd.write(pack_i_0)  # 4 NAMD writes NSTEP or ISTART - NSAVC here?
            dcd.write(pack_ix4_0x4)  # 5, 6, 7, 8
            dcd.write(pack("f", timestep))  # 9 timestep
            dcd.write(pack("i", int(self._unitcell)))  # 10 with unitcell
            dcd.write(pack_ix4_0x4)  # 11, 12, 13, 14
            dcd.write(pack_ix4_0x4)  # 15, 16, 17, 18
            dcd.write(pack("i", 24))  # 19 Pretend to be CHARMM version 24
            dcd.write(pack_i_84)
            dcd.write(pack_i_164)
            dcd.write(pack_i_2)
            dcd.write(b"Created by ProDy".ljust(80))
            temp = now().strftime("%d %B, %Y at %H:%M")
            try:
                temp = bytes(temp, encoding="utf-8")
            except TypeError:
                pass
            dcd.write((b"REMARKS Created " + temp).ljust(80))
            dcd.write(pack_i_164)

            dcd.write(pack_i_4)
            dcd.write(pack(b"i", n_atoms))
            dcd.write(pack_i_4)
            self._first_byte = dcd.tell()

        if self._unitcell:
            if unitcell is None:
                raise TypeError("unitcell data is expected")
            else:
                uc = unitcell
                uc[3:] = np.sin((PISQUARE / 90) * (90 - uc[3:]))
                # uc = uc[[0,2,5,1,3,4]] # <- gpta ??
                uc = uc[[0, 5, 1, 4, 3, 2]]
                pack_i_48 = pack("i", 48)
        dcd.seek(0, 2)

        # for xyz in coords:
        xyz = coords
        if self._unitcell:
            dcd.write(pack_i_48)
            uc.tofile(dcd)
            dcd.write(pack_i_48)
        xyz = xyz.T
        dcd.write(pack_i_4N)
        xyz[0].tofile(dcd)
        dcd.write(pack_i_4N)
        dcd.write(pack_i_4N)
        xyz[1].tofile(dcd)
        dcd.write(pack_i_4N)
        dcd.write(pack_i_4N)
        xyz[2].tofile(dcd)
        dcd.write(pack_i_4N)
        self._n_csets += 1
        dcd.seek(8, 0)
        dcd.write(pack("i", self._n_csets))
        dcd.seek(0, 2)

    def readCoord(self):
        """Read the header information from a dcd file.
        Input: fd - a file struct opened for binary reading.
        Output: 0 on success, negative error code on failure.
        Side effects: *natoms set to number of atoms per frame
                      *nsets set to number of frames in dcd file
                      *istart set to starting timestep of dcd file
                      *nsavc set to timesteps between dcd saves
                      *delta set to value of trajectory timestep
                      *nfixed set to number of fixed atoms
                      *freeind may be set to heap-allocated space
                      *reverse set to one if reverse-endian, zero if not.
                      *charmm set to internal code for handling charmm data.
        """
        LOGGER = logging.getLogger("dynamicEntropy")

        dcd = self._file
        endian = b""  #'=' # native endian
        rec_scale = RECSCALE32BIT
        charmm = None
        dcdcordmagic = unpack(endian + b"i", b"CORD")[0]
        # Check magic number in file header and determine byte order
        bits = dcd.read(calcsize("ii"))

        try:
            temp = unpack(endian + b"ii", bits)
        except:
            return None

        if temp[0] + temp[1] == 84:
            LOGGER.debug("Detected CHARMM -i8 64-bit DCD file of native endianness.")
            rec_scale = RECSCALE64BIT
        elif temp[0] == 84 and temp[1] == dcdcordmagic:
            pass
            LOGGER.debug("Detected standard 32-bit DCD file of native endianness.")
        else:
            if unpack(b">ii", bits) == temp:
                endian = ">"
            else:
                endian = "<"
            temp = unpack(endian + b"ii", bits)
            if temp[0] + temp[1] == 84:
                rec_scale = RECSCALE64BIT
                LOGGER.debug(
                    "Detected CHARMM -i8 64-bit DCD file of opposite endianness."
                )
            else:
                endian = ""
                temp = unpack(endian + b"ii", bits)
                if temp[0] == 84 and temp[1] == dcdcordmagic:
                    LOGGER.debug(
                        "Detected standard 32-bit DCD file of opposite endianness."
                    )
                else:
                    raise IOError("Unrecognized DCD header or unsupported DCD format.")

        # check for magic string, in case of long record markers
        if rec_scale == RECSCALE64BIT:
            raise IOError("CHARMM 64-bit DCD files are not yet supported.")
            temp = unpack(b"I", dcd.read(calcsize("I")))
            if temp[0] != dcdcordmagic:
                raise IOError(
                    "Failed to find CORD magic in CHARMM -i8 64-bit DCD file."
                )

        # Buffer the entire header for random access
        bits = dcd.read(80)

        # CHARMm-genereate DCD files set the last integer in the
        # header, which is unused by X-PLOR, to its version number.
        # Checking if this is nonzero tells us this is a CHARMm file
        # and to look for other CHARMm flags.
        temp = unpack(endian + b"i" * 20, bits)

        if temp[-1] != 0:
            charmm = True

        if charmm:
            LOGGER.debug("CHARMM format DCD file (also NAMD 2.1 and later).")
            temp = unpack(endian + b"i" * 9 + b"f" + b"i" * 10, bits)
        else:
            LOGGER.debug(
                "X-PLOR format DCD file (also NAMD 2.0 and earlier) is not supported."
            )
            return None

        # Store the number of sets of coordinates (NSET)
        self._n_csets = temp[0]
        # Store ISTART, the starting timestep
        self._first_ts = temp[1]
        # Store NSAVC, the number of timesteps between dcd saves
        self._framefreq = temp[2]
        # Store NAMNF, the number of fixed atoms
        self._n_fixed = temp[8]

        if self._n_fixed > 0:
            raise IOError("DCD files with fixed atoms is not yet supported.")

        # Read in the timestep, DELTA
        # Note: DELTA is stored as double with X-PLOR but as float with CHARMm
        self._timestep = temp[9]
        self._unitcell = temp[10] == 1

        # Get the end size of the first block
        if unpack(endian + b"i", dcd.read(rec_scale * calcsize("i")))[0] != 84:
            raise IOError("Unrecognized DCD format.")

        # Read in the size of the next block
        temp = unpack(endian + b"i", dcd.read(rec_scale * calcsize("i")))

        if (temp[0] - 4) % 80 != 0:
            raise IOError("Unrecognized DCD format.")
        noremarks = temp[0] == 84

        # Read NTITLE, the number of 80 character title strings there are
        temp = unpack(endian + b"i", dcd.read(rec_scale * calcsize("i")))

        self._dcdtitle = dcd.read(80)

        if not noremarks:
            self._remarks = dcd.read(80)

        # Get the ending size for this block
        temp = unpack(endian + b"i", dcd.read(rec_scale * calcsize("i")))

        if (temp[0] - 4) % 80 != 0:
            raise IOError("Unrecognized DCD format.")

        # Read in an integer '4'
        if unpack(endian + b"i", dcd.read(rec_scale * calcsize("i")))[0] != 4:
            raise IOError("Unrecognized DCD format.")

        # Read in the number of atoms
        self._n_atoms = unpack(endian + b"i", dcd.read(rec_scale * calcsize("i")))[0]
        # Read in an integer '4'
        if unpack(endian + b"i", dcd.read(rec_scale * calcsize("i")))[0] != 4:
            raise IOError("Bad DCD format.")

        self._is64bit = rec_scale == RECSCALE64BIT
        self._endian = endian
        self._n_floats = (self._n_atoms + 2) * 3

        if self._is64bit:
            if self._unitcell:
                self._bytes_per_frame = 56 + self._n_floats * 8
            else:
                self._bytes_per_frame = self._n_floats * 8
            LOGGER.warning(
                "Reading of 64 bit DCD files has not been tested. "
                "Please report any problems that you may find."
            )
            self._dtype = np.float64
            self._itemsize = 8
        else:
            if self._unitcell:
                self._bytes_per_frame = 56 + self._n_floats * 4
            else:
                self._bytes_per_frame = self._n_floats * 4
            self._dtype = np.float32
            self._itemsize = 4

        self._first_byte = self._file.tell()
        n_csets = (getsize(self._filename) - self._first_byte) // self._bytes_per_frame
        if n_csets != self._n_csets:
            LOGGER.warning(
                "DCD header claims {0} frames, file size "
                "indicates there are actually {1} frames.".format(
                    self._n_csets, n_csets
                )
            )
            self._n_csets = n_csets

        frames = []
        try:
            from alive_progress import alive_it

            bar = alive_it(range(self._n_csets))
        except:
            bar = [x for x in range(self._n_csets)]

        for i in bar:
            unitcell = self._nextUnitcell()
            coords = self._nextCoordset()
            frames.append(Frame(self._nfi, self._n_atoms, coords, unitcell))
        return frames

    def _nextCoordset(self):
        n_floats = self._n_floats
        n_atoms = self._n_atoms
        xyz = fromstring(self._file.read(self._itemsize * n_floats), self._dtype)
        if len(xyz) != n_floats:
            return None
        xyz = xyz.reshape((3, n_atoms + 2)).T[1:-1, :]
        xyz = xyz.reshape((n_atoms, 3))
        self._nfi += 1
        if self._astype is not None and self._astype != xyz.dtype:
            xyz = xyz.astype(self._astype)

        return xyz

    def _nextUnitcell(self):
        if self._unitcell:
            self._file.read(4)
            unitcell = fromstring(self._file.read(48), dtype=np.float64)
            unitcell = unitcell[[0, 2, 5, 1, 3, 4]]
            if np.all(abs(unitcell[3:]) <= 1):
                # This file was generated by CHARMM, or by NAMD > 2.5, with the angle */
                # cosines of the periodic cell angles written to the DCD file.        */
                # This formulation improves rounding behavior for orthogonal cells    */
                # so that the angles end up at precisely 90 degrees, unlike acos().   */
                unitcell[3:] = 90.0 - np.arcsin(unitcell[3:]) * 90 / PISQUARE
            self._file.read(4)
            return unitcell
