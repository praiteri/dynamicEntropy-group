import logging
import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import openmm_wrapper as my


def scaleCharges(system, listOfAtoms, scale=None):
    if scale is None:
        raise Exception("No scaling factor for charges was provided")

    factor = float(scale)
    for f in system.getForces():
        # Rigid ion forcefield
        if isinstance(f, mm.NonbondedForce):
            for idx in listOfAtoms:
                qq = f.getParticleParameters(idx)[0].value_in_unit(
                    unit.elementary_charge
                )
                f.setParticleParameters(idx, qq * factor, 1.0, 0.0)

            # Scale the charge-charge 1-4 interactions
            for ix in range(0, f.getNumExceptions()):
                parm = f.getExceptionParameters(ix)
                if (parm[0] in listOfAtoms) and (parm[1] in listOfAtoms):
                    parm[2] *= factor**2
                elif (parm[0] in listOfAtoms) or (parm[1] in listOfAtoms):
                    parm[2] *= factor
                f.setExceptionParameters(ix, *parm)

        # AMOEBA
        elif isinstance(f, mm.AmoebaMultipoleForce):
            for idx in listOfAtoms:
                parm = f.getMultipoleParameters(idx)
                for i in [0, 1, 2, 9]:
                    parm[i] *= factor
                f.setMultipoleParameters(idx, *parm)
    return


def changeMasses(system, listOfAtoms, listOfMasses):
    if isinstance(listOfMasses, str):
        listOfMasses = [float(x) for x in listOfMasses.split(",")]
    elif isinstance(listOfMasses, float):
        listOfMasses = [listOfMasses]
    assert len(listOfAtoms) == len(
        listOfMasses
    ), "Numer of seleted atoms and number of new masses are different"

    for i, m in zip(listOfAtoms, listOfMasses):
        system.setParticleMass(i, m)

    return


def changeAtomsProperties(system, topology, config):
    logger = logging.getLogger("dynamicEntropy")
    logger.critical("#--- Free energy perturbation -------------#")
    my.debugDictionary(config, name="Modify FF")

    if "charges" in config:
        listOfAtoms = my.getAtoms(config["charges"], topology)
        scaleCharges(system, listOfAtoms, config["charges"]["scale"])

    if "masses" in config:
        listOfAtoms = my.getAtoms(config["masses"], topology)
        changeMasses(system, listOfAtoms, config["masses"]["weight"])
