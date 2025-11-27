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

import logging

import numpy as np
import pandas as pd
import os
import copy


def forwardFluxSampling(setup, modeller, system):

    collectiveVariablesForces = my.createCollectiveVariables(
        setup.config["forwardFluxSampling"]["cv"], modeller.topology
    )
    ffs = FFS(setup, system, collectiveVariablesForces)  # Do not move

    integrator = my.integrator(system, setup.config["md"])
    setup.dumpParametersMD()

    # simulation = app.Simulation(modeller.topology, system, integrator,
    #                             mm.Platform.getPlatformByName(setup.config["basic"]["Platform"]),
    #                             setup.properties)
    simulation = my.createSimulation(
        modeller.topology,
        system,
        integrator,
        mm.Platform.getPlatformByName(setup.config["basic"]["Platform"]),
        setup.properties,
    )

    if setup.config["forwardFluxSampling"]["lambdaID"] == 1:
        simulation.context.setPositions(modeller.positions)
        simulation.context.setVelocitiesToTemperature(
            setup.config["md"]["temperature"],
            setup.rng.integers(low=1, high=99999, size=1),
        )
        simulation.reporters.append(app.StateDataReporter(**setup.output["logOutput"]))

        ffs.computeInitialFlux(simulation, collectiveVariablesForces)
    else:
        ffs.computeCrossingProbability(simulation, collectiveVariablesForces)
    return


class FFS(object):
    def __init__(self, setup, system, cvForce):
        self._logger = logging.getLogger("dynamicEntropy")

        self._rng = np.random.default_rng()

        self.config = setup.config["forwardFluxSampling"]

        self._natoms = system.getNumParticles()
        self._ndof = my.getNumberOfDegreesOfFreedom(system)
        self._scaleFactor = np.sqrt(self._ndof / 3 / (self._natoms - 1))
        self._kBT = setup.config_kBT.value_in_unit(unit.kilojoule_per_mole)

        # Langevin
        m = np.array(
            [
                system.getParticleMass(i).value_in_unit(unit.dalton)
                for i in range(system.getNumParticles())
            ]
        )
        self._masses = np.tile(m[:, np.newaxis], 3)

        dt = setup.config["md"]["timestep"]
        self._sigma = self.config["sigma"] / unit.picoseconds
        self._Aparm = np.exp(-self._sigma * dt)
        self._Bparm = (
            np.sqrt((1 - (self._Aparm) ** 2) * self._kBT / self._masses)
            * self._scaleFactor
        )

        self._cv = mm.CustomCVForce("CV")
        self._cv.setName("FFS_force")
        self._cv.addCollectiveVariable("CV", cvForce)

        system.addForce(self._cv)

        self._lambdaValues = self.config["lambda"]

        if self._lambdaValues[1] > self._lambdaValues[0]:
            self.reachedNextInterface = self.reachedNextFwd
            self.backToReactants = self.backToReactantsFwd
        else:
            self.reachedNextInterface = self.reachedNextBwd
            self.backToReactants = self.backToReactantsBwd

        self._lambdaID = self.config["lambdaID"]
        assert self._lambdaID < len(self._lambdaValues), "Already at the last interface"

        self._lambdaDir = "./lambda_%d" % (self._lambdaID)
        if self._lambdaID > 1:
            self._lambdaDir_prev = "./lambda_%d" % (self._lambdaID - 1)

        if not os.path.isdir(self._lambdaDir):
            os.mkdir(self._lambdaDir)

        self._numberOfSteps = self.config["nsteps"]
        self._samplingFrequency = self.config["nsample"]

        try:
            self._numberOfTrialRuns = self.config["ntrials"]
        except:
            self._numberOfTrialRuns = None

        self._outputFiles = self.fileHandles(self._lambdaDir, self._lambdaID)

        self._logger.info("  kBT = {}".format(self._kBT))
        self._logger.info("  collective variable: {}".format(self.config["cv"]))
        self._logger.info("  lambdaID = {}".format(self._lambdaID))
        self._logger.info("  working directory = {}".format(self._lambdaDir))
        if self._lambdaID > 1:
            self._logger.info("  firing directory = {}".format(self._lambdaDir_prev))
        self._logger.info("  initial lambda = {}".format(self._lambdaValues[0]))
        try:
            self._logger.info(
                "  previous lambda = {}".format(self._lambdaValues[self._lambdaID - 1])
            )
        except:
            pass
        self._logger.info(
            "  next lambda = {}".format(self._lambdaValues[self._lambdaID])
        )
        self._logger.info("  number of steps = {}".format(self._numberOfSteps))
        self._logger.info("  sampling frequency = {}".format(self._samplingFrequency))
        if self._lambdaID > 1:
            self._logger.info(
                "  number of trial runs = {}".format(self._numberOfTrialRuns)
            )

        return

    def __del__(self):
        self._outputFiles["out"].close()
        del self._outputFiles["pos"]
        del self._outputFiles["vel"]
        self._outputFiles["erg"].close()
        try:
            self._outputFiles["rnd"].close()
        except:
            pass

    def getTemperature(self, K):
        kB = (unit.MOLAR_GAS_CONSTANT_R).value_in_unit(
            unit.kilojoule_per_mole / unit.kelvin
        )
        temp = 2 * K / self._ndof / kB
        return temp

    def reachedNextFwd(self, a, b):
        return a >= b

    def backToReactantsFwd(self, a, b):
        return a < b

    def reachedNextBwd(self, a, b):
        return a <= b

    def backToReactantsBwd(self, a, b):
        return a > b

    def fileHandles(self, lambdaDir, lambdaID):
        runID = np.random.randint(0x7FFFFFFF)
        fout = open(os.path.join(lambdaDir, "ffs_%d.out" % (runID)), "a")
        fpos = my.dcdWriter(os.path.join(lambdaDir, "pos_%d.dcd" % (runID)), "wb")
        fvel = my.dcdWriter(os.path.join(lambdaDir, "vel_%d.dcd" % (runID)), "wb")
        ferg = open(os.path.join(lambdaDir, "erg_%d.csv" % (runID)), "a")
        fdbg = open(os.path.join(lambdaDir, "debug_%d.csv" % (runID)), "a")
        ferg.write("U,K,E,T,T0,T1\n")

        if lambdaID > 1:
            frnd = open(os.path.join(lambdaDir, "rnd_%d.out" % (runID)), "a")
        else:
            frnd = None
        handles = {
            "out": fout,
            "pos": fpos,
            "vel": fvel,
            "erg": ferg,
            "rnd": frnd,
            "dbg": fdbg,
        }
        return handles

    def dumpLog0(self, fout, cv, n, t):
        x = float(n) / t
        string = "Initial flux = {:e} (Hz) {} [{}/{}]".format(x, cv, n, t)
        fout.write("%g %g %d ( %g Hz)\n" % (t, cv, n, x))
        fout.flush()
        return string

    def dumpLog(self, fout, nNext, nFire, value, idx):
        string = "Crossing percentage = {:.1f} [{} | {} | {}]".format(
            100 * nNext / nFire, nNext, nFire, idx
        )
        fout.write(
            "%d / %d : %g [%d] ( %g )\n" % (nNext, nFire, value, idx, nNext / nFire)
        )
        fout.flush()
        return string

    def dumpCoordinatesDCD(self, fcoord, fveloc, context):
        state = context.getState(
            getPositions=True, getVelocities=True, enforcePeriodicBox=True
        )
        pos = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)

        cell = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.angstrom)
        unitcell = np.zeros(6, dtype=float)
        unitcell[0] = cell[0][0]
        unitcell[1] = cell[1][1]
        unitcell[2] = cell[2][2]
        unitcell[3] = (
            np.arccos(np.dot(cell[1], cell[2]) / unitcell[1] / unitcell[2])
            * 180.0
            / np.pi
        )
        unitcell[4] = (
            np.arccos(np.dot(cell[0], cell[2]) / unitcell[0] / unitcell[2])
            * 180.0
            / np.pi
        )
        unitcell[5] = (
            np.arccos(np.dot(cell[0], cell[1]) / unitcell[0] / unitcell[1])
            * 180.0
            / np.pi
        )
        fcoord.writeCoord(pos, unitcell)

        vel = state.getVelocities(asNumpy=True)
        fveloc.writeCoord(vel)

    def getDataFromInterface(self, lambdaDir):
        self._logger.info("Reading initial configurations ...")

        ids = [
            f[:-4].split("_")[-1]
            for f in u.findFilesMatchingPattern(lambdaDir, "pos*dcd")
        ]

        startingCoordinates = []
        for idx in ids:
            filename = "pos_" + idx + ".dcd"
            self._logger.info("  {}".format(filename))
            dcd = my.dcdWriter(os.path.join(lambdaDir, filename), "rb")
            newCoord = dcd.readCoord()
            if newCoord is not None:
                startingCoordinates = startingCoordinates + newCoord

        startingVelocities = []
        for idx in ids:
            filename = "vel_" + idx + ".dcd"
            self._logger.info("  {}".format(filename))
            dcd = my.dcdWriter(os.path.join(lambdaDir, filename), "rb")
            newVel = dcd.readCoord()
            if newVel is not None:
                startingVelocities = startingVelocities + newVel

        energiesDF = pd.DataFrame()
        for idx in ids:
            filename = "erg_" + idx + ".csv"
            erg = pd.read_csv(os.path.join(lambdaDir, filename))
            energiesDF = pd.concat([energiesDF, erg], ignore_index=True)

        return startingCoordinates, startingVelocities, energiesDF

    def computeInitialFlux(self, simulation, cv):
        self._logger.info("Computing initial flux ...")

        n = int(self._numberOfSteps / self._samplingFrequency)
        try:
            from alive_progress import alive_it

            bar = alive_it(range(int(n)))
        except:
            bar = [x for x in range(int(n))]

        fcv = open("cv.out", "w")
        goingForward = True
        successfulAttempts = 0

        for i in bar:
            simulation.step(self._samplingFrequency)
            value = self._cv.getCollectiveVariableValues(simulation.context)[0]

            fcv.write("%g\n" % (value))
            fcv.flush()
            if goingForward and self.reachedNextInterface(
                value, self._lambdaValues[self._lambdaID]
            ):
                goingForward = False
                successfulAttempts += 1
                t = simulation.context.getTime().value_in_unit(unit.second)
                string = self.dumpLog0(
                    self._outputFiles["out"], value, successfulAttempts, t
                )
                self._logger.info(string)
                self.dumpCoordinatesDCD(
                    self._outputFiles["pos"],
                    self._outputFiles["vel"],
                    simulation.context,
                )

                s = simulation.context.getState(getEnergy=True)
                U = s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
                K = s.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
                T = self.getTemperature(K)
                self._outputFiles["erg"].write("{},{},{}\n".format(U, K, U + K, T))
                self._outputFiles["erg"].flush()

            elif self.backToReactants(value, self._lambdaValues[0]):
                goingForward = True

        t = simulation.context.getTime().value_in_unit(unit.second)
        self.dumpLog0(self._outputFiles["out"], 0.0, successfulAttempts, t)
        self._logger.info(
            "Initial flux = {:e} (Hz)".format(float(successfulAttempts) / t)
        )

        return

    def computeCrossingProbability(self, simulation, cv):

        self._startingCoordinates, self._startingVelocities, self._energiesDF = (
            self.getDataFromInterface(self._lambdaDir_prev)
        )

        self._logger.info("Computing transition probability ...")
        numFrames = len(self._startingCoordinates)
        self._logger.info(
            "  {:40s} = {}".format("number of initial configurations", numFrames)
        )
        assert numFrames > 0, "No configurations for the firing runs"
        assert len(self._startingCoordinates) == len(
            self._startingVelocities
        ), "Number of frames in coord and vel files are different"

        n = self._numberOfTrialRuns
        try:
            from alive_progress import alive_it

            bar = alive_it(range(int(n)))
        except:
            bar = [x for x in range(int(n))]

        # Loop over number of desired firing attempts
        iRandom = self._rng.integers(0, numFrames, size=self._numberOfTrialRuns)

        np.savetxt(self._outputFiles["rnd"], iRandom, fmt="%d")

        nNext = 0
        for iFire in bar:
            # self.test(simulation)
            idx = iRandom[iFire]
            res, value = self.firingRun(idx, simulation)
            nNext += res
            if res > 0:
                string = self.dumpLog(
                    self._outputFiles["out"], nNext, iFire + 1, value, idx
                )
                self._logger.info(string)

        self._logger.info(
            "Final crossing percentage = {} [{}/{}]".format(
                100 * nNext / (iFire + 1), nNext, iFire + 1
            )
        )

        return

    def firingRun(self, idx, simulation):
        # get energy (potential and kinetics)
        ergInterface = self._energiesDF.iloc[idx, :]

        # get cell and positions
        cell = copy.copy(self._startingCoordinates[idx]._unitcell)
        # import openmm.app as app
        # help(app.internal.unitcell.computeLengthsAndAngles)
        # help(app.internal.unitcell.computePeriodicBoxVectors)
        # help(app.internal.unitcell.reducePeriodicBoxVectors)
        hmat = app.internal.unitcell.computePeriodicBoxVectors(*cell)
        posInterface = copy.copy(self._startingCoordinates[idx]._coords) / 10

        # get velocities
        velInterface = copy.copy(self._startingVelocities[idx]._coords)

        # add cell and positions to the context
        simulation.context.setPeriodicBoxVectors(*hmat)
        simulation.context.setPositions(posInterface)

        # generate noise for the velocities - like a langevin termostat
        shape = self._startingCoordinates[0]._coords.shape
        noise = self._rng.normal(0, 1, size=shape)
        velNew = self._Aparm * velInterface + self._Bparm * noise
        simulation.context.setVelocities(velNew)

        # get new kinetic energy
        s = simulation.context.getState(getEnergy=True)
        U = s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        K = s.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)

        assert (
            np.abs(U - ergInterface["U"]) / ergInterface["U"] < 1e-5
        ), "The potential energy computed for frame {} differs from the old one ({}/{})".format(
            idx, U, ergInterface["U"]
        )

        T0 = self.getTemperature(K)
        T1 = self.getTemperature(ergInterface["K"])
        self._outputFiles["dbg"].write(
            "{},{},{},{}\n".format(K, ergInterface["K"], T0, T1)
        )
        self._outputFiles["dbg"].flush()

        simulation.context.setStepCount(0)
        simulation.context.setTime(0)

        result = 0
        for iStep in range(0, int(self._numberOfSteps / self._samplingFrequency)):
            simulation.step(self._samplingFrequency)
            value = self._cv.getCollectiveVariableValues(simulation.context)[0]

            if self.reachedNextInterface(value, self._lambdaValues[self._lambdaID]):
                # t = simulation.context.getTime().value_in_unit(unit.picosecond)
                self.dumpCoordinatesDCD(
                    self._outputFiles["pos"],
                    self._outputFiles["vel"],
                    simulation.context,
                )

                s = simulation.context.getState(getEnergy=True)
                U = s.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
                K = s.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole)
                T = self.getTemperature(K)
                self._outputFiles["erg"].write(
                    "{},{},{},{},{},{}\n".format(U, K, U + K, T, T0, T1)
                )
                self._outputFiles["erg"].flush()

                result = 1
                break

            elif self.backToReactants(value, self._lambdaValues[0]):
                break
        return result, value

    def test(self, simulation):
        numFrames = len(self._startingCoordinates)
        for idx in range(5):  # numFrames):
            posNew = copy.copy(self._startingCoordinates[idx]._coords) / 10
            velNew = copy.copy(self._startingVelocities[idx]._coords)
            simulation.context.setPositions(posNew)
            simulation.context.setVelocities(velNew)

            s = simulation.context.getState(getEnergy=True)
            U = s.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            K = s.getKineticEnergy().value_in_unit(unit.kilojoules_per_mole)
            # app.PDBFile.writeFile(simulation.topology, 10*posNew, open("test.pdb", 'w'))

        quit()
