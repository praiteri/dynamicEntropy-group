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


def singlePointEnergy(setup, modeller, system):
    # logger = logging.getLogger("dynamicEntropy")

    # simulation = app.Simulation(modeller.topology, system, mm.VerletIntegrator(0),
    #                         mm.Platform.getPlatformByName(setup.config["basic"]["Platform"]),
    #                         setup.properties)

    simulation = my.createSimulation(
        modeller.topology,
        system,
        mm.VerletIntegrator(0),
        mm.Platform.getPlatformByName(setup.config["basic"]["Platform"]),
        setup.properties,
    )

    simulation.context.setPositions(modeller.positions)

    my.computeEnergyFromContext(simulation.context)

    # if setup.debug:
    #     my.computeEnergyFromContext(simulation.context,outputUnit=unit.md_ev)

    return


def rerunTrajectory(setup, modeller, system):
    logger = logging.getLogger("dynamicEntropy")

    my.addElectronVoltToUnit()

    energyUnits = unit.md_ev
    coordUnits = unit.angstrom
    forceUnits = energyUnits / coordUnits

    out1 = open("erg.dat", "w")
    out2 = open("pos.dat", "w")
    out3 = open("force.dat", "w")

    # simulation = app.Simulation(modeller.topology, system, mm.VerletIntegrator(0),
    #                         mm.Platform.getPlatformByName(setup.config["basic"]["Platform"]),
    #                         setup.properties)
    simulation = my.createSimulation(
        modeller.topology,
        system,
        mm.VerletIntegrator(0),
        mm.Platform.getPlatformByName(setup.config["basic"]["Platform"]),
        setup.properties,
    )

    # pdb, cell = my.read_pdb_coordinates('coord.pdb')
    # numberOfAtoms = len(pdb)
    # coords = [x["coord"] for x in pdb]

    # frames = [u.Frame(0,numberOfAtoms,coords,unitcell=cell)]
    numberOfAtoms = modeller.topology.getNumAtoms()
    frames = []
    numberOfFrames = 0

    # Read other frames from DCD
    logger.info("RERUN: Reading the DCD ...")
    dcd = u.dcdWriter("trajectory.0.dcd", "rb")
    pos = dcd.readCoord()

    # Append to frames
    for p in pos:
        frames.append(
            u.Frame(numberOfFrames, numberOfAtoms, p._coords, unitcell=p._unitcell)
        )
        numberOfFrames += 1

    nAtoms = 1440
    out1.write(f"# Energy units = {energyUnits}\n")
    out2.write(f"{nAtoms}\n\n")
    out3.write(f"{nAtoms}\n\n")

    logger.info("RERUN: Processing the frames ...")
    try:
        from alive_progress import alive_it

        bar = alive_it(range(int(numberOfFrames)))
    except:
        bar = [x for x in range(int(numberOfFrames))]

    for f in bar:
        import numpy as np

        pos = np.array(frames[f].getPositions()) * unit.angstrom

        simulation.context.setPositions(pos)
        state = simulation.context.getState(
            getPositions=True, getForces=True, getEnergy=True
        )

        e = state.getPotentialEnergy().value_in_unit(energyUnits)
        forces = state.getForces(asNumpy=True)
        forces = forces.value_in_unit(forceUnits)
        positions = state.getPositions(asNumpy=True).value_in_unit(coordUnits)

        out1.write("%g\n" % (e))
        out1.flush()

        for i in range(nAtoms):
            if i % 5 == 0:
                idx = "Ca"
            elif i % 5 == 1:
                idx = "C4"
            else:
                idx = "O4"

            out2.write(
                "%s %g %g %g\n"
                % (idx, positions[i][0], positions[i][1], positions[i][2])
            )
            out3.write(
                "%s %g %g %g\n" % (idx, forces[i][0], forces[i][1], forces[i][2])
            )

        out2.flush()
        out3.flush()

    return
