### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import openmm as mm
from openmm import CustomIntegrator
from openmm.unit import MOLAR_GAS_CONSTANT_R
import openmm.unit as unit
import new_openmm_wrapper as my

import math

import logging


class CSVRIntegrator(CustomIntegrator):
    """
    CVSR intergrator written by Blake I. Armstrong - Curtin University
    Ref:
        Bussi, G., Donadio, D. & Parrinello, M.
        Canonical Sampling Through Velocity Rescaling.
        Journal of Chemical Physics 126, 014101 (2007).
    """

    def __init__(self, temperature, taut, timestep, system, **kwargs):
        logger = logging.getLogger("dynamicEntropy")
        CustomIntegrator.__init__(self, timestep)

        # Number of degrees of freedom in the system
        self.ndof = my.get_number_of_degrees_of_freedom(system)
        self.nn = self.ndof - 1
        self.kT = MOLAR_GAS_CONSTANT_R * temperature

        self.factor = math.exp(-timestep.value_in_unit(unit.picosecond) / taut)
        self.max_d = kwargs.get("max_d", None)
        if self.max_d is not None:
            try:
                self.max_d = float(self.max_d)
            except ValueError:
                raise ValueError(
                    f"max_d should be a float or convertible to float. Received {self.max_d} of type {type(self.max_d).__name__}"
                )

        logger.debug("------------------------------------#")
        logger.debug("CSVR: kT {}".format(self.kT))
        logger.debug("CSVR: Scale factor {}".format(self.factor))
        logger.debug("------------------------------------#")

        # Integrator initialization.
        self.addGlobalVariable("kT", self.kT)  # thermal energy
        self.addGlobalVariable("ndof", self.ndof)  # number of degrees of freedom
        self.addGlobalVariable("nn", 0)
        self.addGlobalVariable("kk", 0)  # tot kinetic energy
        self.addGlobalVariable("kk2", 0)  # new kinetic energy
        self.addGlobalVariable("sigma", 0)
        self.addGlobalVariable("rr", 0)
        self.addGlobalVariable("rr2", 0)
        self.addGlobalVariable("done", 0)  # gamdev()
        self.addGlobalVariable("v1", 0)  # gamdev()
        self.addGlobalVariable("v2", 0)  # gamdev()
        self.addGlobalVariable("am", 0)  # gamdev()
        self.addGlobalVariable("s", 0)  # gamdev()
        self.addGlobalVariable("a", 0)  # gamdev()
        self.addGlobalVariable("b", 0)  # gamdev()
        self.addGlobalVariable("e", 0)  # gamdev()
        self.addGlobalVariable("check", 0)  # gamdev()
        self.addGlobalVariable("scale", 1.0)
        self.addPerDofVariable("x1", 0)  # for constraints
        if self.max_d is not None:
            self.addPerDofVariable("dd", 0)  # for max_displ

        if self.nn == 0:
            self.addComputeGlobal("nn", "0")
        elif self.nn == 1:
            self.addComputeGlobal("rr2", "gaussian")
            self.addComputeGlobal("nn", "rr2*rr2")
        elif self.nn % 2 == 0:
            self.ia = 0.5 * self.nn
            self.addGlobalVariable("ia", self.ia)
        else:
            self.ia = 0.5 * (self.nn - 1)
            self.addComputeGlobal("rr2", "gaussian")
            self.addGlobalVariable("ia", self.ia)

        self.addComputeSum("kk", "0.5*m*v*v")  # tot kinetic energy
        self.addComputeGlobal("sigma", "ndof*kT*0.5")

        ## gamdev() function
        if self.nn < 1:
            raise ValueError
        if self.nn < 6:
            self.addComputeGlobal("a", "1.0")
            for i in range(self.ia):
                self.addComputeGlobal("a", "a*uniform")
            self.addComputeGlobal("a", "-log(a)")
        else:
            self.addComputeGlobal("done", "0")
            self.beginWhileBlock("done < 1")
            self.addComputeGlobal("v1", "uniform")
            self.addComputeGlobal("v2", "2.0*uniform-1.0")
            self.addComputeGlobal("am", "ia - 1")
            self.addComputeGlobal("s", "sqrt(2.0*am+1.0)")
            self.addComputeGlobal("b", "v2/v1")
            self.addComputeGlobal("a", "s*b+am")
            self.addComputeGlobal("e", "(1+b^2)*exp(am*log(a/am)-s*b)")
            self.addComputeGlobal("check", "v1*v1 + v2*v2")
            self.beginIfBlock("check <= 1.0")
            self.beginIfBlock("a > 0.0")
            self.beginIfBlock("uniform <= e")
            self.addComputeGlobal("done", "done + 1")
            self.endBlock()
            self.endBlock()
            self.endBlock()
            self.endBlock()
        ####

        if self.nn % 2 == 0 and self.nn > 1:
            self.addComputeGlobal("nn", "2.0*a")
        elif self.nn % 2 == 1 and self.nn > 1:
            self.addComputeGlobal("nn", "2.0*a + rr2*rr2")

        self.addGlobalVariable("factor", self.factor)
        self.addComputeGlobal("rr", "gaussian")
        self.addComputeGlobal(
            "kk2",
            "kk+(1.0-factor)*(sigma*(nn+rr*rr)/ndof-kk) + 2.0*rr*sqrt(kk*sigma/ndof*(1.0-factor)*factor)",
        )

        ##re-scale the velocities based on new KE
        self.addComputeGlobal("scale", "sqrt(kk2/kk)")
        self.addComputePerDof("v", "scale*v")

        # Velocity Verlet step
        self.addUpdateContextState()
        self.addComputePerDof("v", "v+0.5*dt*f/m")
        if self.max_d is None:
            self.addComputePerDof("x", "x+dt*v")
        if self.max_d is not None:
            self.addComputePerDof("dd", "dt*v")
            self.addComputePerDof("x", f"x+min({self.max_d},abs(dd))*(2*step(dd)-1)")
        self.addComputePerDof("x1", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt")
        self.addConstrainVelocities()
