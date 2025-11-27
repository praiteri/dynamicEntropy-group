import sys
import openmm as mm
import openmm.app as app
import openmm.unit as unit

# 1. Read PDB
pdb = app.PDBFile("coord.pdb")

# 2. Read forcefield file
forcefield = app.ForceField("oFF.xml")

# 3. Create the system object
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=app.PME,  # Particle Mesh Ewald for long-range electrostatics
    nonbondedCutoff=0.9 * unit.nanometer,  # Cutoff for non-bonded interactions
    constraints=None,  # Constrain bonds involving hydrogen
    rigidWater=False,  # Use rigid water model
    ewaldErrorTolerance=1e-5,  # Error tolerance for PME
    switchDistance=None,
)

# 4. Create force groups for debugging
for n, f in enumerate(system.getForces()):
    f.setForceGroup(n + 1)

# 5. Choose the platform for the simultion
platform = mm.Platform.getPlatformByName("OpenCL")  # Use 'CPU', 'CUDA', or 'OpenCL'

# 6. Create integrator
timestep = 0.001 * unit.picosecond
trelax = 1.0 / unit.picosecond
temperature = 300 * unit.kelvin
pressure = 1 * unit.bar

# integrator = mm.VerletIntegrator(timestep)
integrator = mm.LangevinMiddleIntegrator(temperature, trelax, timestep)
# system.addForce(mm.MonteCarloBarostat(pressure, temperature, 25))

# 7. Create the simulation object
properties = None
simulation = app.Simulation(pdb.topology, system, integrator, platform, properties)

# 8. Add positions and initialise velocities
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temperature, 192837)

# 9. Print energies
print("-" * 5, "Printing Energy Terms", "-" * 5)
state = simulation.context.getState(getPositions=True, getEnergy=True)
c = simulation.context
for force in system.getForces():
    e = c.getState(getEnergy=True, groups={force.getForceGroup()}).getPotentialEnergy()
    print(force.getName(), e)
e = c.getState(getEnergy=True).getPotentialEnergy()
print("Total Energy", e)

# 10. Run MD
print("-" * 5, "Running MD", "-" * 5)
nsteps = 1000
simulation.reporters.append(
    app.StateDataReporter(
        file=sys.stdout,
        reportInterval=100,
        separator="\t",
        step=False,
        time=True,
        potentialEnergy=True,
        kineticEnergy=False,
        totalEnergy=False,
        temperature=True,
        volume=True,
        density=False,
        progress=True,
        remainingTime=True,
        speed=True,
        totalSteps=nsteps,
        elapsedTime=True,
    )
)
simulation.step(nsteps)
