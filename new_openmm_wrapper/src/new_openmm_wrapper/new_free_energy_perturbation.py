### Copyright (C) 2023  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import sys
import copy

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import new_openmm_wrapper as my


class FreeEnergyPerturbationReporter(object):
    """
    Reporter class for Free Energy Perturbation calculations.

    This reporter evaluates the system energy at multiple lambda values during
    the simulation and writes the results to file for post-processing with
    methods like MBAR or thermodynamic integration.
    """

    def __init__(
        self,
        reportInterval,
        context,
        runID=None,
        calculation_type={},
        lambda0=None,
        aux_context=None,
        temperature=300,
    ):
        """
        Initialize the FEP reporter.

        Args:
            reportInterval: Number of steps between energy evaluations
            context: OpenMM simulation context
            filename: Output file for FEP energies (default: 'fep.out')
            lambda0: Current simulation lambda value
            lambda_values: List of lambda values to evaluate energies at
            aux_context: Optional auxiliary context for energy evaluation
        """
        self._reportInterval = reportInterval
        self.context = context
        self.aux_context = aux_context

        self.nvalues = 0
        self.sum_exp = [0, 0]
        self.kBT = (unit.MOLAR_GAS_CONSTANT_R).value_in_unit(
            unit.kilojoules_per_mole / unit.kelvin
        ) * temperature._value

        # Validate required parameters
        if len(calculation_type) == 0:
            my.pretty_log("Simulation type must be provided", logger="error")
            sys.exit(1)
        self.calculation_type = calculation_type

        if runID is None:
            self.output_file_bar = "bar.out"
            self.output_file_mbar = "mbar.out"
            self.output_file_fep = "fep.out"
        else:
            self.output_file_bar = f"bar.{runID}.out"
            self.output_file_mbar = f"mbar.{runID}.out"
            self.output_file_fep = f"fep.{runID}.out"

        if lambda0 is None:
            my.pretty_log("Simulation lambda must be provided", logger="error")
            sys.exit(1)
        self.lambda0 = lambda0

        my.pretty_log(title="FEP reporter", sep=True, logger="debug")
        my.pretty_log({"Report interval": self._reportInterval}, logger="debug")
        my.pretty_log({"kBT (kJ/mol)": self.kBT}, logger="debug")
        my.pretty_log({"λ0": self.lambda0}, logger="debug")
        my.pretty_log(
            {
                "FEP output(s)": [
                    self.output_file_bar,
                    self.output_file_mbar,
                    self.output_file_fep,
                ]
            },
            logger="debug",
        )

    def __del__(self):
        """Close output file when reporter is destroyed."""
        if hasattr(self, "_out_bar"):
            self._out_bar.close()
        if hasattr(self, "_out_mbar"):
            self._out_mbar.close()
        if hasattr(self, "_out_fep"):
            self._out_fep.close()
        pass

    def describeNextReport(self, simulation):
        """
        Describe the next report this reporter will generate.

        Returns:
            Tuple of (steps_until_report, need_positions, need_velocities,
                     need_forces, need_energies, extra_info)
        """
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, False, False, None)

    def report(self, simulation, state):
        """
        Generate a report by evaluating energies at all lambda values.

        Args:
            simulation: The OpenMM Simulation object
            state: Current simulation state
        """
        # Open output file on first report
        for k, v in self.calculation_type.items():
            if k == "mbar" and not hasattr(self, "_out_mbar"):
                str_header = ["Time"] + [f"E{x:5.3f}" for x in v]
                self._out_mbar = open(self.output_file_mbar, "w")
                self._out_mbar.write(" ".join(str_header) + "\n")

            if k == "bar" and not hasattr(self, "_out_bar"):
                str_header = "Time E0 Ef Eb dEf dEb dGf dGb".split(" ")
                self._out_bar = open(self.output_file_bar, "w")
                self._out_bar.write(" ".join(str_header) + "\n")

            if k == "custom" and not hasattr(self, "_out_fep"):
                str_header = [f"Time E{self.lambda0}"] + [f"E{x:5.3f}" for x in v]
                self._out_fep = open(self.output_file_fep, "w")
                self._out_fep.write(" ".join(str_header) + "\n")

        stime = state.getTime().value_in_unit(unit.picosecond)

        if self.aux_context is None:
            for calc, lambda_values in self.calculation_type.items():
                # Evaluate energy at each lambda value
                if calc == "mbar":
                    fep_energies = []
                    for ll in lambda_values:
                        self.context.setParameter("lambda", ll)
                        fep_energies.append(
                            my.computeEnergyFromContext(self.context, quiet=True)
                        )
                    # Always restore original lambda value
                    self.context.setParameter("lambda", self.lambda0)
                    # Write output
                    fep_output = f"{stime:10.3f} " + " ".join(
                        f"{e:15.5f}" for e in fep_energies
                    )
                    self._out_mbar.write(fep_output + "\n")
                    self._out_mbar.flush()

                if calc == "bar":
                    # fep_energies = [my.computeEnergyFromContext(self.context, quiet=True)]
                    fep_energies = []
                    for ll in lambda_values:
                        self.context.setParameter("lambda", ll)
                        fep_energies.append(
                            my.computeEnergyFromContext(self.context, quiet=True)
                        )

                    # Update statistics
                    self.nvalues += 1
                    dE_f = fep_energies[1] - fep_energies[0]
                    dE_b = fep_energies[2] - fep_energies[0]
                    self.sum_exp[0] += np.exp(-dE_f / self.kBT)
                    self.sum_exp[1] += np.exp(-dE_b / self.kBT)

                    # Write BAR output
                    dG_f = -self.kBT * np.log(self.sum_exp[0] / self.nvalues)
                    dG_b = -self.kBT * np.log(self.sum_exp[1] / self.nvalues)

                    fep_output = (
                        f"{stime:10.3f} {fep_energies[0]:15.5f} "
                        f"{fep_energies[1]:15.5f} {fep_energies[2]:15.5f} "
                        f"{dE_f:15.5f} {dE_b:15.5f} {dG_f:15.5f} {dG_b:15.5f}"
                    )
                    # Always restore original lambda value
                    self.context.setParameter("lambda", self.lambda0)
                    # Write output
                    self._out_bar.write(fep_output + "\n")
                    self._out_bar.flush()

                if calc == "custom":
                    fep_energies = [
                        my.computeEnergyFromContext(self.context, quiet=True)
                    ]
                    for ll in lambda_values:
                        self.context.setParameter("lambda", ll)
                        fep_energies.append(
                            my.computeEnergyFromContext(self.context, quiet=True)
                        )
                    # Always restore original lambda value
                    self.context.setParameter("lambda", self.lambda0)
                    # Write output
                    fep_output = f"{stime:10.3f} " + " ".join(
                        f"{e:15.5f}" for e in fep_energies
                    )
                    self._out_fep.write(fep_output + "\n")
                    self._out_fep.flush()


class FreeEnergyPerturbationConstructor(object):
    """
    Main class for constructing and running Free Energy Perturbation simulations.

    Supports multiple FEP types:
    - 'coul': Coulombic/electrostatic interactions
    - 'vdw': Van der Waals interactions
    - 'swap': Particle swap with charge transfer
    - 'swap_qt': Particle swap with charge transfer
    - 'custom': Custom force field modifications
    """

    def __init__(self, config):
        """
        Initialize FEP simulation from configuration dictionary.

        Args:
            config: Configuration dictionary containing simulation parameters
        """
        self.config = config
        self.setup = my.simulationSetup(config)

        # Lambda parameters
        self.variable_name = "lambda"
        self.lambda0 = config["fep"].get("lambda_sim", 0.0)

        self._get_calculation_type()

        # FEP calculation type
        self.fep_type = self.setup.config["fep"].get("type").lower()

        # Create system and integrator
        self.modeller, self.system = my.createSystem(self.setup)
        self.integrator = my.createIntegrator(self.system, self.setup.config["md"])

        # Define atom sets for FEP
        self.all_atoms = set(range(self.modeller.topology.getNumAtoms()))

        if config["fep"].get("type") in ["swap", "swap_qt"]:
            self.fep_atoms1 = set(
                my.get_atoms_from_topology(
                    {"select": config["fep"]["select1"]}, self.modeller.topology
                )
            )
            self.fep_atoms2 = set(
                my.get_atoms_from_topology(
                    {"select": config["fep"]["select2"]}, self.modeller.topology
                )
            )
            self.fep_atoms = self.fep_atoms1 | self.fep_atoms2
        elif config["fep"].get("type") == "custom":
            self.fep_atoms = set()
        else:
            self.fep_atoms = set(
                my.get_atoms_from_topology(config["fep"], self.modeller.topology)
            )

        self.non_fep_atoms = self.all_atoms - self.fep_atoms

        # Calculate equilibration steps
        self._calculate_equilibration_steps()

    def _calculate_equilibration_steps(self):
        """
        Calculate number of equilibration steps from config.

        Can be specified either as:
        - equilibrationSteps: direct number of steps
        - equilibrationTime: time duration (converted to steps using timestep)

        If equilibrationTime is provided, it overrides equilibrationSteps.
        """
        self.eqsteps = self.config["fep"].get("equilibrationSteps", 0)

        eqt = self.config["fep"].get("equilibrationTime", None)
        if eqt is not None:
            # Convert equilibration time to appropriate units
            if isinstance(eqt, str):
                eqt = float(eqt) * unit.picosecond
            elif isinstance(eqt, (int, float)):
                eqt *= unit.picosecond

            # Get timestep from MD config
            if not isinstance(self.setup.config["md"]["timestep"], unit.Quantity):
                dt = self.setup.config["md"]["timestep"] * unit.femtosecond
            else:
                dt = self.setup.config["md"]["timestep"]

            # Calculate steps: equilibrationTime / timestep
            self.eqsteps = int(eqt / dt)

    def _get_calculation_type(self):
        self.calculation_type = {}
        my_dict = self.config["fep"].get("lambda_eval", None)
        if my_dict.get("mbar"):
            dlambda = my_dict.get("mbar")
            if isinstance(dlambda, str):
                dlambda = eval(dlambda)
            elif not isinstance(dlambda, float):
                raise Exception("dlambda must be a float or an evaluatable string")
            lvalues = list(round(x, 4) for x in np.arange(0, 1 + dlambda / 2, dlambda))
            if not any([np.isclose(self.lambda0, x) for x in lvalues]):
                print(lvalues)
                my.pretty_log(
                    lvalues,
                    title=f"λ0 ({self.lambda0}) is not in the list of λ values:",
                    logger="error",
                )
                sys.exit(1)
            self.calculation_type["mbar"] = lvalues

        if my_dict.get("bar"):
            dlambda = my_dict.get("bar")
            lvalues = [
                self.lambda0,
                min(1, self.lambda0 + dlambda),
                max(0, self.lambda0 - dlambda),
            ]
            self.calculation_type["bar"] = lvalues

        if my_dict.get("values"):
            values = my_dict.get("values")
            if isinstance(values, str):
                lvalues = [float(x) for x in values.split(",")]
            elif not isinstance(values, list):
                raise Exception(
                    "lambda values must be a list or a comma separated list"
                )
            else:
                lvalues = [float(x) for x in values]
            self.calculation_type["custom"] = lvalues

    def dump_fep_info(self):
        """Print FEP configuration information for debugging."""
        my.pretty_log(
            self.config["fep"], logger="DEBUG", title="FEP input commands:", sep=True
        )
        my.pretty_log(title="Free Energy Perturbation calculation", sep=True)
        my.pretty_log({"FEP type": self.fep_type}, indent=1)
        my.pretty_log(list(self.fep_atoms), title="FEP atoms:", indent=1, ncols=5)
        my.pretty_log({"Simulation lambda": self.lambda0}, indent=1)
        for k, v in self.calculation_type.items():
            my.pretty_log({f"λ values for {k}": v}, indent=1)

        parameters = self.context.getParameters()
        if len(parameters) > 0:
            my.pretty_log(
                dict(parameters),
                title="Global parameters in context:",
                indent=1,
                logger="DEBUG",
            )

    def create_simulation(self):
        """
        Create the appropriate FEP simulation based on the specified type.

        Dispatches to specific creation methods based on fep_type.
        """
        fep_type_map = {
            "coul": self.create_fep_q_simulation,
            "coulomb": self.create_fep_q_simulation,
            "vdw": self.create_fep_vdw_simulation,
            "vdwl": self.create_fep_vdw_simulation,
            "van_der_waals": self.create_fep_vdw_simulation,
            "custom": self.create_basic_lambda_simulation,
            "swap": self.create_fep_swap_simulation,
            "swap_qt": self.create_fep_swap_simulation,
        }

        fep_creation_method = fep_type_map.get(self.fep_type)
        if fep_creation_method is None:
            my.pretty_log(f"Unknown FEP type: {self.fep_type}", logger="error")
            sys.exit(1)

        fep_creation_method()

        self.simulation = my.createSimulation(
            self.modeller.topology,
            self.system,
            self.integrator,
            mm.Platform.getPlatformByName(self.setup.config["basic"]["Platform"]),
            self.setup.properties,
        )

        self.context = self.simulation.context

    def initialise_md(self):
        """Initialize molecular dynamics simulation with reporters and initial conditions."""
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
        """
        Create a simulation with custom lambda-dependent forces.

        Searches for CustomExternalForce objects with lambda parameters
        and sets them to the specified lambda0 value.
        """
        my.pretty_log(title="Creating custom FEP simulation", sep=True)

        for force in self.system.getForces():
            if isinstance(force, mm.CustomExternalForce):
                # Search for lambda global parameters
                for param_id in range(force.getNumGlobalParameters()):
                    param_name = force.getGlobalParameterName(param_id)
                    if "lambda" in param_name:
                        self.variable_name = param_name
                        force.setGlobalParameterDefaultValue(param_id, self.lambda0)

    def create_fep_q_simulation(self):
        """
        Create a FEP simulation for Coulombic/electrostatic interactions.

        Modifies NonbondedForce to make charges lambda-dependent using
        particle parameter offsets. At lambda=0, charges are zero; at
        lambda=1, charges have their full values.
        """
        my.pretty_log(title="Creating FEP Coulomb simulation", sep=True)

        for force in self.system.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.addGlobalParameter(self.variable_name, self.lambda0)

                for index in self.fep_atoms:
                    charge, sigma, epsilon = force.getParticleParameters(index)
                    my.pretty_log(
                        f"Modifying particle {index}: q={charge}, σ={sigma}, ε={epsilon}",
                        indent=1,
                        logger="DEBUG",
                    )

                    # Set base charge to zero and add lambda-dependent offset
                    # Effective charge = 0 + lambda * charge
                    force.setParticleParameters(index, 0, sigma, epsilon)
                    force.addParticleParameterOffset(
                        self.variable_name, index, charge, 0, 0
                    )

    def create_fep_vdw_simulation(self):
        """
        Create a VDW (Van der Waals) FEP simulation.

        Creates separate force groups for FEP and non-FEP VDW interactions.
        Charges are set to zero for FEP atoms to avoid double-counting
        electrostatics.
        """
        my.pretty_log(title="Creating VDW FEP simulation", sep=True)

        forces_to_remove = []

        for index, force in enumerate(self.system.getForces()):
            # Remove charges from FEP atoms in nonbonded force
            if isinstance(force, mm.NonbondedForce):
                for atom_idx in self.fep_atoms:
                    my.pretty_log(
                        f"Removing charge from atom {atom_idx}", logger="debug"
                    )
                    charge, sigma, epsilon = force.getParticleParameters(atom_idx)
                    force.setParticleParameters(atom_idx, 0, sigma, epsilon)

            # Handle custom nonbonded forces
            if isinstance(force, mm.CustomNonbondedForce):
                # Create non-FEP interactions
                self.add_copy_of_custom_nbf_nofep(
                    copy.deepcopy(force), [self.fep_atoms, self.non_fep_atoms]
                )
                # Create FEP interactions
                self.add_copy_of_custom_nbf_fep(
                    copy.deepcopy(force),
                    "lambda",
                    "(lambda)",
                    [self.fep_atoms, self.non_fep_atoms],
                    force_expr="fep",
                )
                forces_to_remove.append(index)

        # Remove original custom nonbonded forces (must be done after iteration)
        for index in reversed(forces_to_remove):
            self.system.removeForce(index)

        # log_str = {}
        # for n, f in enumerate(self.system.getForces()):
        #     f.setForceGroup(n + 1)
        #     log_str[n + 1] = f.getName()
        # my.pretty_log(title="Creating force groups:", data=log_str, align_width=0, logger="debug")

    def create_fep_swap_simulation(self):
        """
        Create a particle swapping FEP simulation.
        """
        if self.fep_type == "swap_qt":
            my.pretty_log(
                title="Creating FEP simulation for particle swap with charge transfer",
                sep=True,
            )
        else:
            my.pretty_log(title="Creating FEP simulation for particle swap", sep=True)
        group1 = self.fep_atoms1  # Particles disappear
        group2 = self.fep_atoms2  # Particles appear

        # Create exclusions between particle swapping groups
        for force in self.system.getForces():
            if isinstance(force, mm.NonbondedForce):
                charge_group1 = np.sum(
                    [force.getParticleParameters(index)[0]._value for index in group1]
                )
                charge_group2 = np.sum(
                    [force.getParticleParameters(index)[0]._value for index in group2]
                )

                if charge_group1 != charge_group2:
                    my.pretty_log(
                        f"Particle swapping will cause a change in the system charge ({charge_group1} /= {charge_group2})",
                        logger="warning",
                    )
                my.pretty_log(
                    "FEP swap - applying Nonbonded Exceptions", logger="debug"
                )
                for i in group1:
                    for j in group2:
                        force.addException(i, j, 0, 0, 0, replace=True)

            if isinstance(force, mm.CustomNonbondedForce):
                my.pretty_log(
                    "FEP swap - applying CustomNonbonded Exclusions", logger="debug"
                )
                for i in group1:
                    for j in group2:
                        force.addExclusion(i, j)

        forces_to_remove = []

        # Add lambda-dependent charge modifications
        for force_index, force in enumerate(self.system.getForces()):
            if isinstance(force, mm.NonbondedForce):
                force.addGlobalParameter(self.variable_name, self.lambda0)

                for index in self.fep_atoms:
                    charge, sigma, epsilon = force.getParticleParameters(index)
                    my.pretty_log(
                        f"Modifying particle {index}: q={charge}, σ={sigma}, ε={epsilon}",
                        indent=1,
                        logger="DEBUG",
                    )

                    # Apply lambda-dependent charge modifications
                    if index in group1:
                        # Charge disappears: charge → 0
                        force.addParticleParameterOffset(
                            self.variable_name, index, -charge, 0, 0
                        )
                    elif index in group2:
                        force.setParticleParameters(index, 0, sigma, epsilon)
                        # Charge appears: 0 → charge
                        if self.fep_type == "swap":
                            force.addParticleParameterOffset(
                                self.variable_name, index, charge, 0, 0
                            )
                        # Charge appears: 0 → 2*charge
                        elif self.fep_type == "swap_qt":
                            force.addParticleParameterOffset(
                                self.variable_name, index, 2 * charge, 0, 0
                            )
                    else:
                        my.pretty_log("Particle index doesn't exist", logger="error")
                        sys.exit(1)

            # Handle custom nonbonded forces
            if isinstance(force, mm.CustomNonbondedForce):
                # Create non-FEP interactions
                self.add_copy_of_custom_nbf_nofep(
                    copy.deepcopy(force), [group1, group2, self.non_fep_atoms]
                )
                # Create FEP interactions
                self.add_copy_of_custom_nbf_fep(
                    copy.deepcopy(force),
                    "lambda",
                    "(1-lambda)",
                    [group1, self.non_fep_atoms],
                    suffix="1",
                    force_expr="fep2",
                )
                self.add_copy_of_custom_nbf_fep(
                    copy.deepcopy(force),
                    "lambda",
                    "(lambda)",
                    [group2, self.non_fep_atoms],
                    suffix="2",
                    force_expr="fep2",
                )
                forces_to_remove.append(force_index)

        # Remove original custom nonbonded forces (must be done after iteration)
        for index in reversed(forces_to_remove):
            self.system.removeForce(index)

        # log_str = {}
        # for n, f in enumerate(self.system.getForces()):
        #     f.setForceGroup(n + 1)
        #     log_str[n + 1] = f.getName()
        # my.pretty_log(
        #     title="Creating force groups:", data=log_str, align_width=0, logger="debug"
        # )

    def add_copy_of_custom_nbf_nofep(self, force, groups):
        """
        Add a copy of CustomNonbondedForce without FEP modifications.

        Used to handle non-FEP interactions separately in VDW calculations.

        Args:
            force: CustomNonbondedForce to copy
            groups: List of atom groups for interaction
        """
        for group in groups:
            force.addInteractionGroup(group, group)

        force.setName("copy_of_" + force.getName())
        self.system.addForce(force)

    def add_copy_of_custom_nbf_fep(
        self, force, fep_variable, ff_variable, groups, suffix=None, force_expr="fep"
    ):
        """
        Add a lambda-dependent copy of CustomNonbondedForce for FEP.

        Modifies the energy function to include lambda dependence for
        soft-core or other FEP potentials.

        Args:
            force: CustomNonbondedForce to modify
            fep_variable: Name of the lambda parameter
            ff_variable: Variable name in the force field function
            groups: List of atom groups for interaction
        """
        customForceFieldsX, _ = my.get_custom_ff_from_file(self.setup.config)
        if customForceFieldsX is None:
            my.pretty_log("FEP for non-custom ff not implemented", logger="error")
            sys.exit(1)

        # Get the appropriate FEP function from force field definition
        fname = force.getName().replace("_vdwl", "")
        if fname not in customForceFieldsX:
            my.pretty_log("Not sure what to do here...", logger="error")
            sys.exit(1)

        force.addGlobalParameter(fep_variable, self.lambda0)

        # Replace energy function with FEP version
        original_function = [x.strip() for x in force.getEnergyFunction().split(";")]
        fep_function = customForceFieldsX[fname][force_expr].format(var=ff_variable)
        fep_function = [x.strip() for x in fep_function.split(";")]

        my.pretty_log(original_function, title="Original function:", ncols=1, indent=1)
        my.pretty_log(fep_function, title="FEP function:", ncols=1)

        force.setEnergyFunction(";".join(fep_function))
        force.addInteractionGroup(groups[0], groups[1])
        if suffix is not None:
            force.setName(force.getName() + "_" + suffix)
        else:
            force.setName(force.getName() + "_fep")
        self.system.addForce(force)

    def run_test(self):
        """
        Test FEP setup by evaluating energies at all lambda values.

        Useful for verifying the system is set up correctly before
        running production simulations.
        """
        quiet = False
        my.pretty_log(title="Testing FEP context energies", sep=True)
        self.simulation.context.setPositions(self.modeller.positions)

        # Evaluate energy at simulation lambda
        self.context.setParameter(self.variable_name, self.lambda0)
        e0 = my.computeEnergyFromContext(
            self.simulation.context, logger="DEBUG", quiet=quiet
        )

        # Evaluate energies at all lambda values
        self.simulation.context.setPositions(self.modeller.positions)
        for k, v in self.calculation_type.items():
            test_result = {}
            if k == "mbar":
                for ll in v:
                    self.context.setParameter(self.variable_name, ll)
                    e = my.computeEnergyFromContext(
                        self.simulation.context, logger="DEBUG", quiet=quiet
                    )
                    test_result[f"Energy for λ={ll:.3f}"] = (
                        f"E={e:14.5f} ; ΔE={e - e0:15.5f}"
                    )
            if k == "bar":
                for ll in v:
                    self.context.setParameter(self.variable_name, ll)
                    e = my.computeEnergyFromContext(
                        self.simulation.context, logger="DEBUG", quiet=quiet
                    )
                    test_result[f"Energy for λ={ll:.3f}"] = (
                        f"E={e:14.5f} ; ΔE={e - e0:15.5f}"
                    )
            if k == "custom":
                for ll in v:
                    self.context.setParameter(self.variable_name, ll)
                    e = my.computeEnergyFromContext(
                        self.simulation.context, logger="DEBUG", quiet=quiet
                    )
                    test_result[f"Energy for λ={ll:.3f}"] = (
                        f"E={e:14.5f} ; ΔE={e - e0:15.5f}"
                    )

            my.pretty_log(sep=True)
            my.pretty_log(test_result, indent=1)

    def run_equilibration(self):
        """Run equilibration phase at the simulation lambda value."""
        my.pretty_log(title="FEP equilibration phase", sep=True)
        my.pretty_log(f"Equilibration steps = {self.eqsteps} steps", indent=1)
        self.simulation.step(self.eqsteps)

    def run_production(self):
        """
        Run production FEP simulation.

        Adds the FEP reporter to collect energies at multiple lambda values
        during the simulation.
        """
        my.pretty_log(title="FEP production phase", sep=True)
        my.pretty_log(
            f"Production steps = {self.setup.config['md']['numberOfSteps']} steps",
            indent=1,
        )

        # Add FEP reporter
        self.simulation.reporters.append(
            FreeEnergyPerturbationReporter(
                self.config["fep"].get(
                    "reportInterval", self.setup.config["md"]["reportInterval"]
                ),
                self.simulation.context,
                runID=self.config["input"].get("runID"),
                calculation_type=self.calculation_type,
                lambda0=self.lambda0,
                temperature=self.setup.config["md"]["temperature"],
            )
        )

        self.simulation.step(self.setup.config["md"]["numberOfSteps"])


def compute_fep(config):
    """
    Main entry point for Free Energy Perturbation calculations.

    Args:
        config: Configuration dictionary containing all simulation parameters

    The function creates an FEP simulation, optionally runs tests, then
    performs equilibration and production runs.
    """
    fep = FreeEnergyPerturbationConstructor(config)

    fep.create_simulation()

    fep.dump_fep_info()

    # Run test mode if requested
    if config["fep"].get("test", False):
        fep.run_test()
        return

    # Run full simulation
    fep.initialise_md()

    fep.run_equilibration()
    fep.run_production()
