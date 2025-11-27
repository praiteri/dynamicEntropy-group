"""
Osmotic Pressure Calculation Module

This module provides the OsmoticPressureReporter class for calculating and
optionally maintaining constant osmotic pressure in molecular dynamics simulations.

Copyright (C) 2023  Paolo Raiteri
Licensed under GNU General Public License v3 or later
"""

import openmm.unit as unit
import openmm_wrapper as my

import os
import sys
import numpy as np
import logging
from pprint import pformat


logger = logging.getLogger("dynamicEntropy")


class OsmoticPressureReporter(object):
    """Calculate osmotic pressure with optional constant pressure maintenance."""

    # Conversion factor from kJ/mole/nm^3 to bar
    PRESSURE_CONVERSION_FACTOR = 16.605778811

    def __init__(self, runID, resName, config, temperature, topology, context):
        """
        Initialize the osmotic pressure reporter.

        Parameters:
            runID (str): The ID of the current run.
            resName (str): Residue name for restraint.
            config (dict): Configuration settings for osmotic pressure calculation.
            temperature (float): System temperature.
            topology (Topology): OpenMM Topology object.
            context (Context): OpenMM Context object.

        Raises:
            KeyError: If geometry or required config keys are missing.
            Exception: If restraint cannot be found or created.
        """

        self._runID = runID
        self._resName = "_" + resName
        self._context = context
        system = self._context.getSystem()

        # Global parameters for restraint force constant and position
        self._gcmd_parm = ["k" + self._resName, "d0" + self._resName]

        # Initialize restraint
        if config.get("geometry", None) is None:
            self._init_external_restraint(config, system)
        else:
            self.create_restraint_force(config, system, topology)

        logger.critical("Initializing osmotic pressure calculation...")
        logger.debug("Configuration: %s", pformat(config))

        # Find and store force group
        self._find_force_group(system)

        # Setup output file
        self._setup_output_file(config)

        # Setup intervals with validation
        self._setup_intervals(config)

        # Initialize geometry-specific parameters
        self._init_geometry(config, context.getState())

        # Initialize temperature constants
        self._init_temperature_constants(temperature)

        # Setup GCMD if specified
        self._setup_gcmd(config, context)

        logger.info("Osmotic pressure reporter initialized successfully")

    def _init_external_restraint(self, config, system):
        """Initialize using an externally defined restraint."""
        try:
            restraint_force = config["restraint"]["name"]
            self._geometry = config["restraint"]["geometry"].upper()
        except (KeyError, TypeError) as e:
            logger.error("Missing restraint configuration: %s", e)
            sys.exit(1)

        # Verify force exists
        force_found = False
        for f in system.getForces():
            if restraint_force == f.getName():
                force_found = True
                break

        if not force_found:
            logger.error("No restraint force found with name: %s", restraint_force)
            sys.exit(1)

        # Verify global parameters exist
        for param in self._gcmd_parm:
            try:
                self._context.getParameter(param)
            except Exception as e:
                logger.error("Global parameter '%s' not found: %s", param, e)
                sys.exit(1)

        self._force = restraint_force
        logger.info("Using external restraint: %s", self._force)

    def _find_force_group(self, system):
        """Find the force group containing the restraint force."""
        self._forceGroup = set()
        for force in system.getForces():
            if force.getName() == self._force:
                self._forceGroup.add(force.getForceGroup())
                self._restraintForce = force
                break

        if not self._forceGroup:
            logger.error("Restraint force '%s' not found in system", self._force)
            sys.exit(1)

        logger.info("Force name: %s", self._force)
        logger.info("Geometry: %s", self._geometry)

    def _setup_output_file(self, config):
        """Initialize output file with error handling."""
        output_file = config.get("output", f"osmotic.{self._runID}.out")
        try:
            self._out = open(output_file, "w")
            logger.info("Output file: %s", output_file)
        except IOError as e:
            logger.error("Cannot open output file '%s': %s", output_file, e)
            sys.exit(1)

    def _setup_intervals(self, config):
        """Setup and validate computation intervals."""
        self._computeInterval = config.get("computeInterval", 1000)
        self._reportInterval = config.get("reportInterval", self._computeInterval)
        self._sampleLength = config.get("sampleLength", 1000)

        # Validate intervals
        if self._reportInterval < self._computeInterval:
            logger.error(
                "reportInterval (%d) < computeInterval (%d)",
                self._reportInterval,
                self._computeInterval,
            )
            sys.exit(1)

        if self._reportInterval % self._computeInterval != 0:
            logger.error("reportInterval must be multiple of computeInterval")
            sys.exit(1)

        logger.info("Compute interval: %d steps", self._computeInterval)
        logger.info("Report interval: %d steps", self._reportInterval)

    def _init_geometry(self, config, state):
        """Initialize geometry-specific parameters."""
        self._dirs = (0, 1, 2, 1)

        if self._geometry.upper() != "SPHERE":
            try:
                axis_config = config["geometry"][self._geometry.lower()]
                self._dir = axis_config["axis"]
            except (KeyError, TypeError):
                logger.error("Missing axis for %s geometry", self._geometry.upper())
                sys.exit(1)

            # Parse axis direction
            vers = 1
            self._dir = self._dir.replace("+", "")
            if "-" in self._dir:
                vers = -1
                self._dir = self._dir.replace("-", "")

            # Map to direction indices
            axis_upper = self._dir.upper()
            if axis_upper == "X":
                self._dirs = (2, 1, 0, vers)
            elif axis_upper == "Y":
                self._dirs = (0, 2, 1, vers)
            elif axis_upper == "Z":
                self._dirs = (0, 1, 2, vers)
            else:
                logger.error(
                    "Unknown crystallographic direction: %s (%s)",
                    self._geometry.upper(),
                    self._dir,
                )
                sys.exit(1)

            logger.info("Axis direction: %s", self._dir.upper())

        self._box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(
            unit.nanometer
        )

        # Initialize geometry-specific area and volume functions
        geometry = self._geometry.upper()

        if geometry == "SPHERE":
            self._compute_area = (
                lambda r: 4
                * np.pi
                * (r**2 + 2 * r * self._factor_1 + 2 * self._factor_0)
            )
            self._compute_volume = (
                lambda r: 4
                * np.pi
                * (
                    r**3 / 3
                    + (r**2 + self._factor_0) * self._factor_1
                    + r * 2 * self._factor_0
                )
            )
            self._update_gcmd_parameter = (
                lambda gcmd_parameter, mu: gcmd_parameter * mu ** (1 / 3)
            )

        elif geometry == "HARMONIC":
            i, j = self._dirs[0], self._dirs[1]
            self._compute_area = lambda r: self._box[i][i] * self._box[j][j]
            self._compute_volume = (
                lambda r: self._box[i][i] * self._box[j][j] * self._factor_1
            )
            self._update_gcmd_parameter = lambda gcmd_parameter, mu: gcmd_parameter

        elif geometry == "SLAB":
            i, j = self._dirs[0], self._dirs[1]
            self._compute_area = lambda r: 2 * self._box[i][i] * self._box[j][j]
            self._compute_volume = (
                lambda r: self._box[i][i] * self._box[j][j] * 2 * (r + self._factor_1)
            )
            self._update_gcmd_parameter = lambda gcmd_parameter, mu: gcmd_parameter * mu

        elif geometry == "PLANE":
            i, j = self._dirs[0], self._dirs[1]
            self._compute_area = lambda r: self._box[i][i] * self._box[j][j]
            self._compute_volume = lambda r: 0
            self._update_gcmd_parameter = (
                lambda gcmd_parameter, mu: gcmd_parameter + (mu - 1) * 10
            )

        elif geometry == "CYLINDER":
            k = self._dirs[2]
            rtmp = self._box[k][k] * 2 * np.pi
            self._compute_area = lambda r: rtmp * (r + self._factor_1)
            self._compute_volume = lambda r: rtmp * (
                0.5 * r**2 + r * self._factor_1 + self._factor_0
            )
            self._update_gcmd_parameter = (
                lambda gcmd_parameter, mu: gcmd_parameter * mu ** (1 / 2)
            )

    def _init_temperature_constants(self, temperature):
        """Initialize temperature-dependent constants."""
        self._T = temperature
        self._kt = (unit.MOLAR_GAS_CONSTANT_R * self._T).in_units_of(
            unit.kilojoules_per_mole
        )
        self._kappa = (
            self._context.getParameter(self._gcmd_parm[0])
            * unit.kilojoules_per_mole
            / unit.nanometer**2
        )
        self._factor_0 = self._kt / (
            self._context.getParameter(self._gcmd_parm[0]) * unit.kilojoules_per_mole
        )
        self._factor_1 = np.sqrt(np.pi * self._factor_0 / 2)

        self._localCounter = 0
        self._averagePressure = 0.0
        self._osmoticPressure = np.zeros(self._sampleLength)

        logger.debug(
            "Temperature: %.1f K, kBT: %.2f kJ/mol", self._T._value, self._kt._value
        )

    def _setup_gcmd(self, config, context):
        """Setup constant osmotic pressure (GCMD) parameters."""
        self._gcmd = config.get("gcmd", None)

        if self._gcmd is None:
            logger.info("Running at constant volume (no GCMD)")
            self._write_output_header(self._gcmd)
            return

        # Set GCMD defaults
        gcmd_defaults = {
            "pext": 0.0,
            "K": 0.01,
            "tau": 1.0,
            "sample": 100,
        }
        for k, v in self._gcmd.items():
            if k not in gcmd_defaults and k != "restart":
                logger.error("Unknown GCMD parameter: %s", k)
                sys.exit(1)
            if k != "restart":
                gcmd_defaults[k] = v

        self._gcmd_pext = float(gcmd_defaults["pext"])
        self._gcmd_compressibility = 0.01  # Can be modified
        self._gcmd_tau = float(gcmd_defaults["tau"])
        self._gcmd_sample = int(gcmd_defaults["sample"])
        self._gcmd_dt = (
            context.getIntegrator().getStepSize().value_in_unit(unit.picosecond)
        )
        self._gmcd_constant = (
            self._gcmd_compressibility * self._gcmd_dt / self._gcmd_tau
        )

        self._osmoticPressure.fill(self._gcmd_pext)

        logger.info("Constant Osmotic Pressure (GCMD) setup:")
        logger.info("  Target pressure: %.2f bar", self._gcmd_pext)
        logger.info("  Compressibility: %.4f 1/bar", self._gcmd_compressibility)
        logger.info("  Tau: %.1f ps", self._gcmd_tau)
        logger.info("  Timestep: %.4f ps", self._gcmd_dt)

        # Handle restart file
        if "restart" in self._gcmd:
            self._load_gcmd_restart()
        else:
            self._gcmd_press = self._gcmd_pext

        self._write_output_header(self._gcmd)

    def _load_gcmd_restart(self):
        """Load GCMD state from restart file."""
        restart_file = self._gcmd["restart"]

        if not os.path.exists(restart_file):
            logger.error("GCMD restart file not found: %s", restart_file)
            sys.exit(1)

        try:
            with open(restart_file, "r") as f:
                lines = f.readlines()

            if not lines:
                logger.error("GCMD restart file is empty: %s", restart_file)
                sys.exit(1)

            values = lines[-1].split()
            self._localCounter = self._gcmd_sample
            self._context.setParameter(self._gcmd_parm[1], float(values[4]))
            self._gcmd_press = float(values[3])

            logger.info("Restarted GCMD from: %s", restart_file)
            logger.info("  Restraint parameter: %s", values[3])
            logger.info("  Average osmotic pressure: %s bar", values[4])

        except (IndexError, ValueError) as e:
            logger.error("Invalid GCMD restart file format: %s", e)
            sys.exit(1)

    def _write_output_header(self, gcmd=None):
        """Write output file header."""
        header = (
            f"# Osmotic Pressure Calculation\n"
            f"# Number of solute particles: {self._restraintForce.getNumParticles()}\n"
            f"# Geometry: {self._geometry}\n"
            f"# Temperature (K): {self._T._value}\n"
            f"# kBT (kJ/mol): {self._kt._value}\n"
            f"# Force constant (kJ/mol/nm^2): {self._kappa._value}\n"
            f"# Computation interval (steps): {self._computeInterval}\n"
            f"# Report interval (steps): {self._reportInterval}\n"
            f"# Sample length (steps): {self._sampleLength}\n"
        )
        if gcmd is not None:
            header += f"# Target osmotic pressure (bar): {self._gcmd_pext}\n"
            header += f"# Compressibility (1/bar): {self._gcmd_compressibility}\n"
            header += f"# Tau (ps): {self._gcmd_tau}\n"

        header_fields = (
            f"# {'Time (ps)':>20s} "
            f"{'Pi(t) (bar)':>20s} "
            f"{'<Pi> (bar)':>20s} "
            f"{'Pi(gcmd) (bar)':>20s} "
            f"{'Parm (nm)':>20s} "
            f"{'Volume (nm^3)':>20s} "
            f"{'Concentration (M)':>20s}\n"
        )

        self._out.write(header)
        self._out.write(header_fields)
        self._out.flush()

    def __del__(self):
        """Cleanup on deletion."""
        try:
            self._out.close()
        except Exception:
            pass

    def describeNextReport(self, simulation):
        """Describe the next report frequency."""
        steps = self._computeInterval - simulation.currentStep % self._computeInterval
        return (steps, False, False, True, False)

    def report(self, simulation, state):
        """Process a report step."""
        self._localCounter += 1
        stime = state.getTime().value_in_unit(unit.picosecond)

        # Extract forces
        forces = simulation.context.getState(
            getForces=True, groups=self._forceGroup
        ).getForces(asNumpy=True)
        forces = forces.value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
        f = np.sum(np.sqrt(np.einsum("ij,ij->i", forces, forces)))

        gcmd_parameter = self._context.getParameter(self._gcmd_parm[1])
        # box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
        # i, j, k = self._dirs[0:3]

        # # Calculate area and volume based on geometry
        # if self._geometry.upper() == "SPHERE":
        #     area = (
        #         4
        #         * np.pi
        #         * (gcmd_parameter**2 + 2 * gcmd_parameter * self._factor_1 + 2 * self._factor_0)
        #     )
        #     volume = (
        #         4
        #         * np.pi
        #         * (
        #             gcmd_parameter**3 / 3
        #             + (gcmd_parameter**2 + self._factor_0) * self._factor_1
        #             + gcmd_parameter * 2 * self._factor_0
        #         )
        #     )

        # elif self._geometry.upper() == "HARMONIC":
        #     area = box[i][i] * box[j][j]
        #     volume = area * self._factor_1

        # elif self._geometry.upper() == "SLAB":
        #     area = box[i][i] * box[j][j]
        #     volume = area * 2 * (gcmd_parameter + self._factor_1)
        #     area *= 2

        # elif self._geometry.upper() == "PLANE":
        #     area = box[i][i] * box[j][j]
        #     volume = 0

        # elif self._geometry.upper() == "CYLINDER":
        #     rtmp = box[k][k] * 2 * np.pi
        #     area = rtmp * (gcmd_parameter + self._factor_1)
        #     volume = rtmp * (0.5 * gcmd_parameter**2 + gcmd_parameter * self._factor_1 + self._factor_0)

        # else:
        #     logger.error("Unknown geometry: %s", self._geometry)
        #     return

        # print(area, volume)
        area = self._compute_area(gcmd_parameter)
        volume = self._compute_volume(gcmd_parameter)

        # Convert kJ/mole/nm^3 to bar
        opress = f / area * self.PRESSURE_CONVERSION_FACTOR
        self._osmoticPressure[(self._localCounter - 1) % self._sampleLength] = opress

        # Update global average pressure
        self._averagePressure = (
            self._averagePressure * (self._localCounter - 1) + opress
        ) / self._localCounter

        # Calculate window average
        if self._gcmd is None:
            n = min(self._localCounter, self._sampleLength)
            self._gcmd_press = np.mean(self._osmoticPressure[:n])
        else:
            # Apply Berendsen-like barostat
            self._gcmd_press = np.mean(self._osmoticPressure)
            mu = 1.0 - self._dirs[3] * self._gmcd_constant * (
                self._gcmd_pext - self._gcmd_press
            )

            gcmd_parameter = self._update_gcmd_parameter(gcmd_parameter, mu)
            # if self._geometry.upper() == "PLANE":
            #     gcmd_parameter = gcmd_parameter + (mu - 1) * 10
            # elif self._geometry.upper() == "SLAB":
            #     gcmd_parameter = gcmd_parameter * mu
            # elif self._geometry.upper() == "CYLINDER":
            #     gcmd_parameter = gcmd_parameter * mu ** (1 / 2)
            # elif self._geometry.upper() == "SPHERE":
            #     gcmd_parameter = gcmd_parameter * mu ** (1 / 3)

            self._context.setParameter(self._gcmd_parm[1], gcmd_parameter)

        # Write output at report interval
        if simulation.currentStep % self._reportInterval != 0:
            return

        try:
            conc = self._restraintForce.getNumParticles() / volume / 0.6022140857
        except ZeroDivisionError:
            conc = 0.0

        self._out.write(
            "  "
            f"{stime:20.1f} "
            f"{opress:20.5f} "
            f"{self._averagePressure:20.5f} "
            f"{self._gcmd_press:20.5f} "
            f"{gcmd_parameter:20.5f} "
            f"{volume:20.5f} "
            f"{conc:20.5f}\n"
        )
        self._out.flush()

    def create_restraint_force(self, config, system, topology):
        """Create a new restraint force from configuration."""
        logger = logging.getLogger("osmotic_pressure")

        if "sphere" in config["geometry"]:
            self._geometry = "SPHERE"
            cmd = config["geometry"]["sphere"]
            my.check_required_keywords(cmd, ["species", "kappa", "radius", "centre"])

            val = ",".join([str(cmd["kappa"]), str(cmd["radius"]), cmd["centre"]])
            expr = f"0.5*{self._gcmd_parm[0]}*(max(0,d-{self._gcmd_parm[1]}))^2;d=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)"
            cmd = {
                f"osmoticWall{self._resName}": {
                    "global": f"{self._gcmd_parm[0]},{self._gcmd_parm[1]}",
                    "par": f"{self._gcmd_parm[0]},{self._gcmd_parm[1]},x0,y0,z0",
                    "species": cmd["species"],
                    "val": val,
                    "fullexp": expr,
                }
            }

        elif "plane" in config["geometry"]:
            self._geometry = "PLANE"
            cmd = config["geometry"]["plane"]
            my.check_required_keywords(cmd, ["species", "kappa", "pos", "axis"])

            val = ",".join([str(cmd["kappa"]), str(cmd["pos"])])
            p = cmd["axis"].replace("-", "").replace("+", "")

            if "-" in cmd["axis"]:
                expr = f"0.5*{self._gcmd_parm[0]}*min(0,{p}-{self._gcmd_parm[1]})^2"
            else:
                expr = f"0.5*{self._gcmd_parm[0]}*max(0,{p}-{self._gcmd_parm[1]})^2"

            cmd = {
                f"osmoticWall{self._resName}": {
                    "global": f"{self._gcmd_parm[0]},{self._gcmd_parm[1]}",
                    "par": f"{self._gcmd_parm[0]},{self._gcmd_parm[1]}",
                    "species": cmd["species"],
                    "val": val,
                    "fullexp": expr,
                }
            }

        elif "slab" in config["geometry"]:
            self._geometry = "SLAB"
            cmd = config["geometry"]["slab"]
            my.check_required_keywords(
                cmd, ["species", "kappa", "width", "centre", "axis"]
            )

            val = ",".join([str(cmd["kappa"]), str(cmd["width"]), str(cmd["centre"])])
            p = cmd["axis"]
            p0 = p + "0"
            expr = (
                f"0.5*{self._gcmd_parm[0]}*(max(0,d-{self._gcmd_parm[1]}))^2;"
                f"d=abs({p}-{p0})"
            )
            cmd = {
                f"osmoticWall{self._resName}": {
                    "global": f"{self._gcmd_parm[0]},{self._gcmd_parm[1]}",
                    "par": f"{self._gcmd_parm[0]},{self._gcmd_parm[1]},{p0}",
                    "species": cmd["species"],
                    "val": val,
                    "exp": expr,
                }
            }

        elif "cylinder" in config["geometry"]:
            self._geometry = "CYLINDER"
            cmd = config["geometry"]["cylinder"]
            my.check_required_keywords(
                cmd, ["species", "kappa", "radius", "centre", "axis"]
            )

            val = ",".join([str(cmd["kappa"]), str(cmd["radius"]), str(cmd["centre"])])

            p = cmd["axis"]
            if p.upper() == "X":
                expr = f"0.5*{self._gcmd_parm[0]}*(max(0,d-{self._gcmd_parm[1]}))^2;d=sqrt((y-y0)^2+(z-z0)^2)"
                p0 = "y0,z0"
            elif p.upper() == "Y":
                expr = f"0.5*{self._gcmd_parm[0]}*(max(0,d-{self._gcmd_parm[1]}))^2;d=sqrt((x-x0)^2+(z-z0)^2)"
                p0 = "x0,z0"
            elif p.upper() == "Z":
                expr = f"0.5*{self._gcmd_parm[0]}*(max(0,d-{self._gcmd_parm[1]}))^2;d=sqrt((x-x0)^2+(y-y0)^2)"
                p0 = "x0,y0"

            cmd = {
                f"osmoticWall{self._resName}": {
                    "global": f"{self._gcmd_parm[0]},{self._gcmd_parm[1]}",
                    "par": f"{self._gcmd_parm[0]},{self._gcmd_parm[1]},{p0}",
                    "species": cmd["species"],
                    "val": val,
                    "fullexp": expr,
                }
            }
        else:
            raise KeyError(
                f"Unknown geometry for Osmotic Pressure calculation: {self._geometry}"
            )

        logger.info("Geometry: %s", self._geometry.lower())

        # Add restraint, need to reinitialise the context
        my.addRestraints(cmd, system, topology)
        self._context.reinitialize(preserveState=True)

        self._force = f"osmoticWall{self._resName}"
