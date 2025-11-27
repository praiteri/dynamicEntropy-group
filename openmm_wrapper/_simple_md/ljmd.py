#!/usr/bin/env python3

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

import os
import ast
import sys
import numpy as np
import logging

import openmm as mm
import openmm.app as app
import openmm.unit as unit

# Simulation parameters
config = {
    "timestep"                   : 0.002 ,
    "numberOfSteps"              : 100000,
    "temperature"                : 150 * unit.kelvin,
    "thermostatParameter"        : 1.0 / unit.picoseconds,
    "NPT"                        : False,
    "pressure"                   : 1 * unit.bar,
    "barostatUpdate"             : 25,
}

# Forcefield paramters for Ar and Kr
atomTypes = [
    {"name":"Ar","mass":39.95,"sigma":0.3405,"epsilon":0.996015},
    {"name":"Kr","mass":39.95,"sigma":0.3405,"epsilon":0.996015},
    # {"name":"Kr","mass":83.80,"sigma":0.3670,"epsilon":1.388420},
]
atomNames = [x["name"] for x in atomTypes]

# Read the coordinates
pdb  = app.PDBFile("coord.pdb")

# Create the system object
system = mm.System()
system.setDefaultPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())
for a in pdb.topology.atoms():
    idx = atomNames.index(a.name)
    system.addParticle(atomTypes[idx]["mass"] * unit.amu)

# Create the forcefield
M = len(atomTypes)
epsilonAR_r = np.zeros((M,M), dtype="float64")
sigmaAR_r = np.zeros((M,M), dtype="float64")

for i in range(M):
    for j in range(i,M):
        epsilonAR_r[i][j] = np.sqrt(atomTypes[i]["epsilon"]*atomTypes[i]["epsilon"])
        epsilonAR_r[j][i] = epsilonAR_r[i][j]
        sigmaAR_r[i][j] = 0.5*(atomTypes[i]["sigma"]+atomTypes[i]["sigma"])
        sigmaAR_r[j][i] = sigmaAR_r[i][j]

# The Lennard-Jones potential we will create in OpenMM accepts these arrays in list form
epsilonLST_r = (epsilonAR_r).ravel().tolist()
sigmaLST_r   = (sigmaAR_r).ravel().tolist()

# Create the force object
customNonbondedForce = mm.CustomNonbondedForce('4*eps*((sig/d)^12-(sig/d)^6); d=max(0.1,r); eps=epsilon(type1, type2); sig=sigma(type1, type2)')
customNonbondedForce.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
customNonbondedForce.setCutoffDistance(1.0 * unit.nanometer)
customNonbondedForce.addTabulatedFunction('epsilon', mm.Discrete2DFunction(M, M, epsilonLST_r))
customNonbondedForce.addTabulatedFunction('sigma', mm.Discrete2DFunction(M, M, sigmaLST_r))
customNonbondedForce.addPerParticleParameter('type')

for a in pdb.topology.atoms():
    idx = atomNames.index(a.name)
    customNonbondedForce.addParticle([idx])

# Add the force to the system
system.addForce(customNonbondedForce)

for n,f in enumerate(system.getForces()):
    f.setForceGroup(n+1)

# Create the integrator object
integrator = mm.LangevinMiddleIntegrator( 
    config["temperature"],
    config["thermostatParameter"],
    config["timestep"])

# Add the barostat for NPT simulation
if config["NPT"]:
    system.addForce(
        mm.MonteCarloBarostat(
            config["pressure"], 
            config["temperature"], 
            config["barostatUpdate"]))

# Create the simulation object
simulation = app.Simulation(pdb.topology, system, integrator,
                            mm.Platform.getPlatformByName("OpenCL"),
                            #mm.Platform.getPlatformByName("CUDA"),
                            {'Precision' : 'mixed'})

# Add the velocities to the simulation
simulation.context.setPositions(pdb.positions)

# Screen output
simulation.reporters.append(
    app.StateDataReporter(
        sys.stdout, 1000, totalSteps = int(config["numberOfSteps"]), separator= "\t",
        step=False, time=True, potentialEnergy=False, kineticEnergy=False,
        totalEnergy=False, temperature=True, volume=False, density=True,
        progress=True, remainingTime=True, speed=True, elapsedTime=False
    )
)

# File output
simulation.reporters.append(
    app.StateDataReporter(
        "md.log" , 1000, separator= ",",
        step=False, time=True, potentialEnergy=True, kineticEnergy=False,
        totalEnergy=False, temperature=True, volume=True, density=True,
        progress=False, remainingTime=False, speed=False, elapsedTime=False
    )
)

# Trajectory output
simulation.reporters.append(
    app.DCDReporter( "traj.dcd" , 1000 )
)

simulation.reporters.append(
    OsmoticPressureReporter(simulation.context))

# Energy minimisation
e = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print(f"Initial energy {e}")
simulation.minimizeEnergy(tolerance=0.001)
e = simulation.context.getState(getEnergy=True).getPotentialEnergy()
print(f"Energy after minimisation {e}")

#s=simulation.context.getState(getPositions=True)
#p=s.getPositions()
#simulation.context.reinitialize()
#simulation.context.setPositions(p)

# Generate the velocities for the simulation
simulation.context.setVelocitiesToTemperature(
    config["temperature"] , 1267 )

# Run molecular dynamics
simulation.step( int(config["numberOfSteps"]) )
