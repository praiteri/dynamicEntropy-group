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

class OsmoticPressureReporter(object):
    def __init__(self, context, config):
        self._force = "Membrane"
        self._forceGroup = set()
        system = context.getSystem()
        for force in system.getForces():
            if force.getName() == self._force:
                self._forceGroup.add(force.getForceGroup())
                break

        f =  "osmotic.out".format()
        self._out = open(f, 'w')

        self._reportInterval  = config["reportInterval"]
        self._computeInterval  = config["computeInterval"]
        self._sampleLength  = config["sampleLength"]
        geom2surf = {"SPHERE" : 0, "PLANE" : 1, "SLAB" : 2, "HARMONIC" : -1}
        self._geometry = config["geometry"].upper()
        self._gcmd_context = context

        try:
            numSurf = geom2surf[self._geometry]
        except KeyError:
            raise KeyError(f"Unknown geometry for Osmotic Pressure calculation: {self._geometry}")
        self._numSurf = numSurf
        if self._geometry != "SPHERE":
            try:
                self._dir = config["axis"]
            except:
                raise Exception(f"Missing axis for {self._geometry} geometry")
                
            if self._dir.upper() == "X":
                self._dirs = (2,1,0)
            elif self._dir.upper() == "Y":
                self._dirs = (0,2,1)
            elif self._dir.upper() == "Z":
                self._dirs = (0,1,2)
            else:
                raise Exception("Unknown crystallographic direction for {} ({})".format(self._geometry.upper(),self._dir))

        self._localCounter = 0
        self._averagePressure = 0.0

        defaults = {
            'pext'   : 0.0 ,
            'parm'   : "d0" ,
            'K'      : 0.01 ,
            'tau'    : 1.0 ,
            'sample' : 100 ,
        }
        self._gcmd_parm   = defaults["parm"]
        self._K      = float(config["parms"]["k"]) 
        self._kt      = float(config["kt"]) 

        if "gcmd" not in config:
            self._osmoticPressure = np.zeros(self._sampleLength)
            self._out.write('#Time Pi(t) <Pi> Pi(sample)\n')
            self._out.flush()
            self._gcmd = False
            return

        self._gcmd = config["gcmd"]

        for k,v in self._gcmd.items():
            try: defaults[k] = v
            except: raise Exception("Unknown GCMD parameter")
            
        self._gcmd_pext   = float(defaults["pext"])
        self._gcmd_K      = float(defaults["K"]) 
        self._gcmd_tau    = float(defaults["tau"])
        self._gcmd_sample = int(defaults["sample"])
        self._gcmd_dt     = context.getIntegrator().getStepSize().value_in_unit(unit.picosecond)

        self._osmoticPressure = np.full(self._sampleLength, self._gcmd_pext)
        

        if "restart" in self._gcmd:
            if not os.path.exists(self._gcmd["restart"]):
                raise Exception("GCMD: missing restart file ({})".format(self._gcmd["restart"]))
            with open(self._gcmd["restart"],"r") as f:
                lines = f.readlines()
            if lines:
                values = lines[-1].split()
            else:
                raise Exception("GCMD: cannot restart from file ({})".format(self._gcmd["restart"])) 
            
            self._localCounter = self._gcmd_sample
            self._gcmd_context.setParameter(self._gcmd_parm,float(values[4]))
            self._gcmd_press = float(values[3])

        self._out.write(f'#Time Pi(t) <Pi> Pi(gcmd) {defaults["parm"]}\n')
        self._out.flush()

        return

    
    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        steps = self._computeInterval - simulation.currentStep%self._computeInterval
        return (steps, False, False, True, False)

    def report(self, simulation, state):
        self._localCounter += 1
        stime = state.getTime().value_in_unit(unit.picosecond)

        # extract forces
        forces = simulation.context.getState(getForces=True, groups=self._forceGroup).getForces(asNumpy=True)
        forces = forces.value_in_unit(unit.kilojoule_per_mole/unit.nanometer)

        # spherical restraining potential
        if self._numSurf == 0:
            self._radius = self._gcmd_context.getParameter(self._gcmd_parm)
            self._radius += np.sqrt(np.pi * self._kt / 4 / self._K)
            area = 4 * np.pi * self._radius**2
            f = np.sum(np.sqrt(np.einsum('ij,ij->i', forces, forces)))

        # Harmonic restraint
        elif self._numSurf == -1:
            box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
            i,j,k = self._dirs
            area = (box[i][i] * box[j][j])
            f = np.sum(forces[:,k])
            
        # flat restraining potential (flat bottom well or plane)
        else:
            box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
            i,j,k = self._dirs
            area = (box[i][i] * box[j][j]) * self._numSurf
            f = np.sum(np.abs(forces[:,k]))

        # convert 'kJ/mole/nm^3' to 'bar'
        opress = f / area * 16.605388
        self._osmoticPressure[ (self._localCounter-1)%self._sampleLength ] = opress

        # Global average of the pressure
        self._averagePressure =  \
            (self._averagePressure * (self._localCounter-1) + opress) / self._localCounter

        if not self._gcmd:
            # Window average of the pressure
            n = min(self._localCounter , self._sampleLength)
            self._gcmd_press = np.mean(self._osmoticPressure[:n])
            if simulation.currentStep%self._reportInterval == 0:
                self._out.write('{} {} {} {}\n'.format(stime, opress, self._averagePressure, self._gcmd_press))
                self._out.flush()

        # GCMD
        else:
            self._gcmd_press = np.mean(self._osmoticPressure)
            # Berendsen-like barostat
            # \mu = 1 - \frac{\kappa_T \delta t}{\tau_\Pi} (\Pi_0 -\Pi)
            # \kappa = 0.01 (bar^{-1}) (~ isothermal compressibilitty)
            # \delta t = simulation timestep
            # \tau_\Pi = 1.0 ps
            mu = 1.0 - self._gcmd_K*self._gcmd_dt/self._gcmd_tau * (self._gcmd_pext-self._gcmd_press)
            parm = mu * self._gcmd_context.getParameter(self._gcmd_parm)
            self._gcmd_context.setParameter(self._gcmd_parm,parm)
    
            if simulation.currentStep%self._reportInterval == 0:
                self._out.write('{} {} {} {} {}\n'.format(stime, opress, self._averagePressure, self._gcmd_press, parm))
                self._out.flush()

        return
    
# Simulation parameters
config = {
    "timestep"                   : 0.002 ,
    "numberOfSteps"              : 1e8,
    "temperature"                : 150 * unit.kelvin,
    "thermostatParameter"        : 1.0 / unit.picoseconds,
    "NPT"                        : False,
    "pressure"                   : 1 * unit.bar,
    "barostatUpdate"             : 25,
    "screenReport"               : 1000,
    "trajReport"                 : 1000,
    "fileReport"                 : 1000,
    "coordinatesFile"           : "coord.pdb"
}

osmoticConfig = {
   "geometry" : "sphere",
   "computeInterval" : 1000,
   "reportInterval" : 1000,
   "sampleLength" : 1000,
   "gcmd" : { 
        "parm" : "d0",  
        "pext" : 30, 
        "tau" : 1.0, 
    },
   "directions" : ("x","y","z"),
   "kt" : 8.3144640540282e-3 * (config["temperature"] / unit.kelvin),
   "parms" : {
       "k" : 100,
       "d0" : 3.0,
       "x0" : 3.5,
       "y0" : 3.5,
       "z0" : 3.5,
    },
}

# Forcefield paramters for Ar and Kr
atomTypes = [
    {"name":"Ar","mass":39.95,"sigma":0.3405,"epsilon":0.996015},
    {"name":"Kr","mass":39.95,"sigma":0.3405,"epsilon":0.996015},
    # {"name":"Kr","mass":83.80,"sigma":0.3670,"epsilon":1.388420},
]
atomNames = [x["name"] for x in atomTypes]

# Read the coordinates
pdb  = app.PDBFile(config["coordinatesFile"])

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

energy = "k*(max(0,d-d0))^2"
distance = " d = periodicdistance(x,y,z,x0,y0,z0)"
if "x" not in osmoticConfig["directions"]: distance = distance.replace("x0","x")
if "y" not in osmoticConfig["directions"]: distance = distance.replace("y0","y")
if "z" not in osmoticConfig["directions"]: distance = distance.replace("z0","z")
expression = energy + " ; " + distance
force = mm.CustomExternalForce(expression)

listOfSoluteSpecies = ["Kr"]
atomsList = [ atom.index for atom in pdb.topology.atoms() if atom.name in listOfSoluteSpecies]

for direction in osmoticConfig["directions"]:
    direction_start = f"{direction}0"
    try:
        p0 = osmoticConfig["parms"][direction_start]
    except KeyError:
        raise ValueError(f"Need a {direction_start} parameter for direction {direction}")
    force.addGlobalParameter(direction_start, p0)
force.addGlobalParameter("d0", osmoticConfig["parms"]["d0"])
force.addGlobalParameter("k", osmoticConfig["parms"]["k"])

for i in atomsList:
    force.addParticle(i,())

force.setName("Membrane")
system.addForce(force)

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
                            #mm.Platform.getPlatformByName("OpenCL"),
                            #mm.Platform.getPlatformByName("CUDA"),
                            mm.Platform.getPlatformByName("HIP"),
                            {'Precision' : 'mixed'})

# Add the velocities to the simulation
simulation.context.setPositions(pdb.positions)

# Screen output
simulation.reporters.append(
    app.StateDataReporter(
        sys.stdout, config["screenReport"], totalSteps = int(config["numberOfSteps"]), separator= "\t",
        step=False, time=True, potentialEnergy=False, kineticEnergy=False,
        totalEnergy=False, temperature=True, volume=False, density=True,
        progress=True, remainingTime=True, speed=True, elapsedTime=False
    )
)

# File output
simulation.reporters.append(
    app.StateDataReporter(
        "md.log" , config["fileReport"], separator= ",",
        step=False, time=True, potentialEnergy=True, kineticEnergy=False,
        totalEnergy=False, temperature=True, volume=True, density=True,
        progress=False, remainingTime=False, speed=False, elapsedTime=False
    )
)

# Trajectory output
simulation.reporters.append(
    app.DCDReporter( "traj.dcd" , config["trajReport"] )
)

simulation.reporters.append(
    OsmoticPressureReporter(simulation.context, osmoticConfig))

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
r = np.random.randint(1, 99999)
simulation.context.setVelocitiesToTemperature(
    config["temperature"] , r )

# Run molecular dynamics
simulation.step( int(config["numberOfSteps"]) )
