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

import logging
import openmm as mm
import openmm.app as app
import openmm_wrapper as my

import pathlib

# class minimizationReporter(mm.MinimizationReporter):
#
#    # within the class you can declare variables that persist throughout the
#    # minimization
#
#    energies = [] # array to record progress
#
#    # you must override the report method and it must have this signature.
#    def report(self, iteration, x, grad, args):
#        '''
#        the report method is called every iteration of the minimization.
#
#        Args:
#            iteration (int): The index of the current iteration. This refers
#                             to the current call to the L-BFGS optimizer.
#                             Each time the minimizer increases the restraint strength,
#                             the iteration index is reset to 0.
#
#            x (array-like): The current particle positions in flattened order:
#                            the three coordinates of the first particle,
#                            then the three coordinates of the second particle, etc.
#
#            grad (array-like): The current gradient of the objective function
#                               (potential energy plus restraint energy) with
#                               respect to the particle coordinates, in flattened order.
#
#            args (dict): Additional statistics described above about the current state of minimization.
#                         In particular:
#                         “system energy”: the current potential energy of the system
#                         “restraint energy”: the energy of the harmonic restraints
#                         “restraint strength”: the force constant of the restraints (in kJ/mol/nm^2)
#                         “max constraint error”: the maximum relative error in the length of any constraint
#
#        Returns:
#            bool : Specify if minimization should be stopped.
#        '''
#
#        # Within the report method you write the code you want to be executed at
#        # each iteration of the minimization.
#        # In this example we get the current energy, print it to the screen, and save it to an array.
#        logger = logging.getLogger('dynamicEntropy')
#
#        current_energy = args['system energy']
#
#        if iteration % 100 == 0: # only print to screen every 100 iterations for clarity of webpage display
#            logger.debug("MINIMIZER: {:5d} {:8.3f}".format(iteration,current_energy))
#
#        self.energies.append(current_energy)
#
#        # The report method must return a bool specifying if minimization should be stopped.
#        # You can use this functionality for early termination.
#        return False


def minimise(setup, modeller, system):
    logger = logging.getLogger("dynamicEntropy")
    setup.dumpParametersMD()

    # simulation = app.Simulation(modeller.topology, system, mm.VerletIntegrator(0),
    #                             mm.Platform.getPlatformByName(setup.config["basic"]["Platform"]),
    #                             setup.properties)
    simulation = my.createSimulation(
        modeller.topology,
        system,
        mm.VerletIntegrator(0),
        mm.Platform.getPlatformByName(setup.config["basic"]["Platform"]),
        setup.properties,
    )

    simulation.context.setPositions(modeller.positions)

    logger.info("Energy minimisation ...")
    my.computeEnergyFromContext(
        simulation.context, header="#--- Initial energy -----------------------#"
    )

    #    reporter=my.minimizationReporter()
    #
    #    simulation.minimizeEnergy(setup.config["minimise"]["tolerance"],
    #                              setup.config["minimise"]["maxIterations"],
    #                              reporter=reporter)

    # print(setup.config["minimise"])
    simulation.minimizeEnergy(
        setup.config["minimise"]["tolerance"], setup.config["minimise"]["maxIterations"]
    )

    my.computeEnergyFromContext(
        simulation.context, header="#--- Final energy -------------------------#"
    )

    positions = simulation.context.getState(getPositions=True).getPositions()
    simulation.context.reinitialize()
    simulation.context.setPositions(positions)

    app.PDBFile.writeFile(
        simulation.topology, positions, open(setup.config["minimise"]["output"], "w")
    )
    logger.critical("#-------------------------------------#")

    # import matplotlib.pyplot as plt
    # plt.plot(reporter.energies)
    # plt.ylabel("System energy (kJ/mol)")
    # plt.xlabel("Minimization iteration")
    # plt.show()
