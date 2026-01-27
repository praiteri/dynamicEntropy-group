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

import sys
import copy

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import new_openmm_wrapper as my


class FreeEnergyPerturbationReporter(object):
    def __init__(
        self,
        reportInterval,
        context,
        filename=None,
        lambda0=None,
        lamda_values=None,
        aux_context=None,
    ):
        self._reportInterval = reportInterval
        self.context = context
        self.aux_context = aux_context

        self.file = filename or "fep.out"

        self.lambda0 = lambda0
        if self.lambda0 is None:
            my.pretty_log("Simulation lambda must be provided", logger="error")
            sys.exit(1)

        self.lamda_values = lamda_values
        if self.lamda_values is None:
            my.pretty_log(
                "Lambda values for FEP force reporting must be provided", logger="error"
            )
            sys.exit(1)

    def __del__(self):
        if hasattr(self, "_out"):
            self._out.close()

    def describeNextReport(self, simulation):
        """
        # Returns a five element tuple.
        # The first element is the number of steps until the next report.
        # The remaining elements specify whether that report will require
        # positions, velocities, forces, and energies respectively
        """
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, False, False, None)

    def report(self, simulation, state):
        if not hasattr(self, "_out"):
            self._out = open(self.file, "w")

        stime = state.getTime().value_in_unit(unit.picosecond)

        if self.aux_context is None:
            # state = simulation.context.getState(getPositions=True, groups=-1)
            # e0 = my.computeEnergyFromContext(simulation.context, quiet=True)
            fep_energies = []
            for ll in self.lamda_values:
                self.context.setParameter("lambda", ll)
                e = my.computeEnergyFromContext(self.context, quiet=True)
                fep_energies.append(f"{e:15.5f}")
                # fep_energies.append(f"{e - e0:15.5f}")
            self._out.write(f"{stime:10.3f}" + " ".join(fep_energies) + "\n")

            # Reset lambda to original value
            self.context.setParameter("lambda", self.lambda0)


class FreeEnergyPerturbationConstructor(object):
    def __init__(self, config):
        self.config = config
        self.setup = my.simulationSetup(config)

        self.variable_name = "lambda"
        self.lambda0 = config["fep"].get("lambda_sim", 0.0)
        self.lambda_values = config["fep"].get("lambda_eval", None)

        self.fep_type = self.setup.config["fep"].get("type").lower()

        self.modeller, self.system = my.createSystem(self.setup)
        self.integrator = my.createIntegrator(self.system, self.setup.config["md"])

        self.all_atoms = set(range(self.modeller.topology.getNumAtoms()))
        self.fep_atoms = set(
            my.get_atoms_from_topology(config["fep"], self.modeller.topology)
        )
        self.non_fep_atoms = self.all_atoms - self.fep_atoms

        self.eqsteps = config["fep"].get("equilibrationSteps", 0)
        eqt = config["fep"].get("equilibrationTime", None)
        if eqt is not None:
            if isinstance(eqt, str):
                eqt = float(eqt) * unit.picosecond
            elif isinstance(eqt, (int, float)):
                eqt *= unit.picosecond

            if not isinstance(self.setup.config["md"]["timestep"], unit.Quantity):
                dt = self.setup.config["md"]["timestep"] * unit.femtosecond
            else:
                dt = self.setup.config["md"]["timestep"]
            self.eqsteps = int(eqt / dt)

    def dump_fep_info(self):
        my.pretty_log(
            self.config["fep"], logger="DEBUG", title="FEP input commands:", sep=True
        )
        my.pretty_log(title="Free Energy Perturbation calculation", sep=True)
        my.pretty_log({"FEP type": self.fep_type}, indent=1)
        my.pretty_log(list(self.fep_atoms), title="FEP atoms:", indent=1, ncols=5)
        my.pretty_log({"Simulation lambda": self.lambda0}, indent=1)
        my.pretty_log({"Evaluation lambda": self.lambda_values}, indent=1)

        parameters = self.context.getParameters()
        if len(parameters) > 0:
            my.pretty_log(
                dict(parameters),
                title="Global parameters in context:",
                indent=1,
                logger="DEBUG",
            )

    def create_simulation(self):
        if self.fep_type in [
            "coul",
            "electrostatics",
            "electrostatic",
            "charge",
            "coulomb",
        ]:
            self.create_fep_q_simulation()

        elif self.fep_type in ["vdw", "van der waals", "van_der_waals", "vdwl"]:
            self.create_fep_vdw_simulation()

        elif self.fep_type in ["custom"]:
            self.create_basic_lambda_simulation()

        elif self.fep_type in ["qt"]:
            self.create_qtransfer_simulation()

        else:
            my.pretty_log(f"Uknown FEP type {self.fep_type}", logger="error")
            sys.exit(1)

        self.context = self.simulation.context

    def initialise_md(self):
        _ = my.simulationReporters(
            self.simulation,
            runID=self.setup.config["input"]["runID"],
            reportInterval=self.setup.config["md"]["reportInterval"],
            numberOfSteps=self.setup.config["md"].get("numberOfSteps", None),
            configReporters=self.setup.config.get("reporters", None),
        )

        my.initialiseMolecularDynamics(
            self.simulation,
            self.modeller.positions,
            self.setup.config["md"]["temperature"],
        )

    def create_basic_lambda_simulation(self):
        my.pretty_log(title="Creating custom FEP simulation", sep=True)
        for force in self.system.getForces():
            if isinstance(force, mm.CustomExternalForce):
                # print(force.getNumGlobalParameters())
                for id in range(force.getNumGlobalParameters()):
                    # print(
                    #     force.getGlobalParameterName(id),
                    #     force.getGlobalParameterDefaultValue(id),
                    #     force.getNumPerParticleParameters(),
                    # )
                    # for idx in range(force.getNumParticles()):
                    #     print(force.getParticleParameters(idx))
                    if "lambda" in force.getGlobalParameterName(id):
                        self.variable_name = force.getGlobalParameterName(id)
                        force.setGlobalParameterDefaultValue(id, self.lambda0)

        self.simulation = my.createSimulation(
            self.modeller.topology,
            self.system,
            self.integrator,
            mm.Platform.getPlatformByName(self.setup.config["basic"]["Platform"]),
            self.setup.properties,
        )

    def create_qtransfer_simulation(self):
        my.pretty_log(title="Creating FEP charge-transfer simulation", sep=True)

        # Create two group to transfer charge
        group1 = [0, 1, 2, 3, 4]  # Charge appears
        group2 = [5]  # Charge disappears
        group3 = []  # Charge is doubled (eg Na+ to Ca2+)

        # Create appropriate exclusions
        for force in self.system.getForces():
            print(force.getName())
            if isinstance(force, mm.NonbondedForce):
                my.pretty_log(
                    "FEP qt - applying Nonbonded Exceptions",
                    logger="debug",
                )
                for i in group1:
                    for j in group2:
                        force.addException(i, j, 0, 0, 1)

            if isinstance(force, mm.CustomNonbondedForce):
                my.pretty_log(
                    "FEP qt - applying CustomNonbonded Exclusions",
                    logger="debug",
                )
                for i in group1:
                    for j in group2:
                        force.addExclusion(i, j)

        # Add global parameter for FEP
        for force in self.system.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.addGlobalParameter(self.variable_name, self.lambda0)
                for index in self.fep_atoms:
                    # Original parameters
                    charge, sigma, epsilon = force.getParticleParameters(index)
                    my.pretty_log(
                        f"Modifying particle {index}: q={charge}, σ={sigma}, ε={epsilon}",
                        indent=1,
                        logger="DEBUG",
                    )
                    # Add an offset that depends on the global parameter
                    # When charge_scale=1, the particle has its normal charge
                    # When charge_scale=0, the particle has zero charge
                    # atom charge = ParticleParameters.charge + charge_scale * ParticleParameterOffset.charge
                    # There are two equivalent ways to do this:

                    if index in group1:
                        force.setParticleParameters(index, 0, sigma, epsilon)
                        force.addParticleParameterOffset(
                            self.variable_name, index, charge, 0, 0
                        )
                    elif index in group2:
                        force.addParticleParameterOffset(
                            self.variable_name, index, -charge, 0, 0
                        )
                    elif index in group3:
                        force.addParticleParameterOffset(
                            self.variable_name, index, charge, 0, 0
                        )

                    else:
                        my.pretty_log("Particle index doesn't exist", logger="error")
                        sys.exit(1)

        customForceFieldsX, parms = self.get_custom_ff_from_file()
        print(customForceFieldsX)
        print(parms)
        quit()
        self.simulation = my.createSimulation(
            self.modeller.topology,
            self.system,
            self.integrator,
            mm.Platform.getPlatformByName(self.setup.config["basic"]["Platform"]),
            self.setup.properties,
        )

    def create_fep_q_simulation(self):
        """
        Create a FEP simulation for Coulomb interactions
        """
        my.pretty_log(title="Creating FEP Coulomb simulation", sep=True)
        for force in self.system.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.addGlobalParameter(self.variable_name, self.lambda0)
                for index in self.fep_atoms:
                    # Original parameters
                    charge, sigma, epsilon = force.getParticleParameters(index)
                    my.pretty_log(
                        f"Modifying particle {index}: q={charge}, σ={sigma}, ε={epsilon}",
                        indent=1,
                        logger="DEBUG",
                    )
                    # Add an offset that depends on the global parameter
                    # When charge_scale=1, the particle has its normal charge
                    # When charge_scale=0, the particle has zero charge
                    # atom charge = ParticleParameters.charge + charge_scale * ParticleParameterOffset.charge
                    # There are two equivalent ways to do this:
                    force.setParticleParameters(index, 0, sigma, epsilon)
                    force.addParticleParameterOffset(
                        self.variable_name, index, charge, 0, 0
                    )
                    # force.addParticleParameterOffset(name, index, -charge, 0, 0)

        self.simulation = my.createSimulation(
            self.modeller.topology,
            self.system,
            self.integrator,
            mm.Platform.getPlatformByName(self.setup.config["basic"]["Platform"]),
            self.setup.properties,
        )

    def add_copy_of_custom_nbf_nofep(self, force, groups):
        for gg in groups:
            force.addInteractionGroup(gg, gg)

        force.setName("copy_of_" + force.getName())
        n = self.system.getNumForces()
        force.setForceGroup(n + 1)
        self.system.addForce(force)

    def add_copy_of_custom_nbf_fep(self, force, fep_variable, ff_variable, groups):
        customForceFieldsX, _ = my.get_custom_ff_from_file(self.setup.config)
        if customForceFieldsX is None:
            my.pretty_log("FEP for noncustom ff not implemented", logger="error")
            sys.exit(1)

        fname = force.getName().replace("_vdwl", "")
        if fname not in customForceFieldsX:
            my.pretty_log("Not sure what to do here...", logger="error")
            sys.exit(1)

        force.addGlobalParameter(fep_variable, self.lambda0)

        original_function = [x.strip() for x in force.getEnergyFunction().split(";")]
        fep_function = customForceFieldsX[fname]["fep"].format(var=ff_variable)
        fep_function = [x.strip() for x in fep_function.split(";")]

        my.pretty_log(original_function, title="Original function:", ncols=1, indent=1)
        my.pretty_log(fep_function, title="FEP function:", ncols=1)

        force.setEnergyFunction(";".join(fep_function))
        force.addInteractionGroup(groups[0], groups[1])
        n = self.system.getNumForces()
        force.setForceGroup(n + 1)
        self.system.addForce(force)

    def create_fep_vdw_simulation(self):

        my.pretty_log(title="Creating VDW FEP simulation", sep=True)

        # Set charges to zero
        for index, force in enumerate(self.system.getForces()):

            # Ideally there should be a test here to decide
            # whether we should break out of the loop

            if isinstance(force, mm.NonbondedForce):
                for index in self.fep_atoms:
                    my.pretty_log(f"Removing charge to atom {index}", logger="debug")
                    charge, sigma, epsilon = force.getParticleParameters(index)
                    force.setParticleParameters(index, 0, sigma, epsilon)

            if isinstance(force, mm.CustomNonbondedForce):
                # Create a copy of the force for the rest of the system that will not be modified
                self.add_copy_of_custom_nbf_nofep(
                    copy.deepcopy(force), [self.fep_atoms, self.non_fep_atoms]
                )
                self.add_copy_of_custom_nbf_fep(
                    copy.deepcopy(force),
                    "lambda",
                    "(lambda)",
                    [self.fep_atoms, self.non_fep_atoms],
                )
                self.system.removeForce(index)

        self.simulation = my.createSimulation(
            self.modeller.topology,
            self.system,
            self.integrator,
            mm.Platform.getPlatformByName(self.setup.config["basic"]["Platform"]),
            self.setup.properties,
        )

    def run_test(self):
        quiet = False
        my.pretty_log(title="Testing FEP context energies", sep=True)
        self.simulation.context.setPositions(self.modeller.positions)

        self.context.setParameter(self.variable_name, self.lambda0)
        e0 = my.computeEnergyFromContext(
            self.simulation.context, logger="DEBUG", quiet=quiet
        )

        self.simulation.context.setPositions(self.modeller.positions)
        test_result = {}
        for ll in self.lambda_values:
            self.context.setParameter(self.variable_name, ll)
            e = my.computeEnergyFromContext(
                self.simulation.context, logger="DEBUG", quiet=quiet
            )
            test_result[f"Energy for λ={ll:.3f}"] = f"E={e:14.5f} ; ΔE={e - e0:15.5f}"
        my.pretty_log(sep=True)
        my.pretty_log(test_result, indent=1)

    def run_equilibration(self):

        my.pretty_log(title="FEP equilibration phase", sep=True)
        my.pretty_log(f"Equilibration steps =  {self.eqsteps} steps", indent=1)
        self.simulation.step(self.eqsteps)

    def run_production(self):
        my.pretty_log(title="FEP production phase", sep=True)
        my.pretty_log(
            f"Production steps =  {self.setup.config['md']['numberOfSteps']} steps",
            indent=1,
        )

        self.simulation.reporters.append(
            FreeEnergyPerturbationReporter(
                self.config["fep"].get(
                    "reportInterval", self.setup.config["md"]["reportInterval"]
                ),
                self.simulation.context,
                self.config["fep"].get("output", "fep.out"),
                lambda0=self.lambda0,
                lamda_values=self.lambda_values,
            )
        )

        self.simulation.step(self.setup.config["md"]["numberOfSteps"])


def compute_fep(config):
    """Free energy perturbation"""

    fep = FreeEnergyPerturbationConstructor(config)
    fep.create_simulation()

    fep.dump_fep_info()

    if config["fep"].get("test", False):
        fep.run_test()
        return

    fep.initialise_md()

    fep.run_equilibration()

    fep.run_production()
