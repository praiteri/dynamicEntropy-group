import logging
import pathlib

import openmm as mm
import openmm.app as app
import openmm.unit as unit
import openmm_wrapper as my


def createSimulation(topology, system, integrator, platform, properties):
    logger = logging.getLogger("dynamicEntropy")

    simulation = app.Simulation(topology, system, integrator, platform, properties)

    try:
        nonbonded = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][
            0
        ]
    except:
        nonbonded = [
            f for f in system.getForces() if isinstance(f, mm.AmoebaMultipoleForce)
        ][0]

    nbm = nonbonded.getNonbondedMethod()
    nbm_list = {
        0: "NoCutoff",
        1: "CutoffNonPeriodic",
        2: "CutoffPeriodic",
        3: "Ewald",
        4: "PME",
        5: "LJPME",
    }

    logger.critical("Simulation details ...")
    logger.info("  {:40s} = {}".format("Nonbonded method", nbm_list[nbm]))
    logger.info(
        "  {:40s} = {}".format("Global cutoff distance", nonbonded.getCutoffDistance())
    )
    try:
        logger.debug(
            "Use switching distance: {}".format(nonbonded.getUseSwitchingFunction())
        )
        logger.debug("Switching distance: {}".format(nonbonded.getSwitchingDistance()))
    except:
        pass

    if nbm in [4, 5]:
        logger.info(
            "  {:40s} = {}".format(
                "PME parameters",
                nonbonded.getPMEParametersInContext(simulation.context),
            )
        )
    if nbm in [3, 4, 5]:
        logger.info(
            "  {:40s} = {}".format(
                "Ewald tolerance", nonbonded.getEwaldErrorTolerance()
            )
        )
        logger.debug(
            "Reciprocal space force group: {}".format(
                nonbonded.getReciprocalSpaceForceGroup()
            )
        )

    if nbm in [2]:
        logger.info(
            "  {:40s} = {}".format(
                "Reaction field dielectric", nonbonded.getReactionFieldDielectric()
            )
        )

    return simulation
