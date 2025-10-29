import sys
import logging

from ..logger import formatting
from ..orca_tools.get_values import get_values
from ..orca_tools.get_coordinates import get_coordinates
from ..orca_tools.get_geopt import get_geopt
from ..orca_tools.get_vibrational_modes import get_vibrational_modes
from ..orca_tools.get_thermochemistry import get_thermochemistry
from ..orca_tools.get_tddft import get_tddft

from ..orca_tools.get_input_commands import get_input_commands
from ..orca_tools.get_charges import get_charges
from ..orca_tools.get_orbital_energies import get_orbital_energies

from ..orca_tools.utils import list_to_kwargs


def run_summary(orca_file):
    """ """
    logger = logging.getLogger("mylogger")

    # Extract general information about the system
    number_of_atoms = get_values(orca_file, "Number of atoms")
    number_of_atoms.write(first=True)

    # Extract the coordinates
    coordinates = get_coordinates(orca_file)
    coordinates.print_formula()

    # Charge and multiplicity
    charge = get_values(orca_file, "Total Charge           Charge", text="Total charge")
    charge.write(first=True)
    multiplicity = get_values(
        orca_file, "Multiplicity           Mult", text="Multiplicity"
    )
    multiplicity.write(last=True)

    # Number of electrons
    nalpha = get_values(orca_file, "N(Alpha)", text="Number of alpha electrons")
    nalpha.write(first=True)
    nbeta = get_values(orca_file, "N(Beta)", text="Number of beta electrons")
    nbeta.write(first=True)
    number_of_electrons = get_values(orca_file, "NEL", text="Total number of Electrons")
    number_of_electrons.write(first=True)

    # Geometry optimization data
    geometry_optimisation = get_geopt(orca_file)
    geometry_optimisation.write()

    # Check the cluster integrity
    coordinates.check_cluster_integrity()

    # Extract the vibrational modes
    vibrational_modes = get_vibrational_modes(orca_file)
    vibrational_modes.write()

    # Thermochemistry
    thermochemistry = get_thermochemistry(orca_file)
    thermochemistry.write()

    # TD-DFT data
    tddft = get_tddft(orca_file)
    if tddft.result is not None:
        tddft.write(nalpha=nalpha.result[0], nbeta=nbeta.result[0])

    logger.info(formatting().dashes)
    logger.info("Done!")


def process_orca_file(args):
    """ """
    mf = formatting()
    logger = logging.getLogger("mylogger")

    # Open the file
    filename = args["input"]
    if filename is None:
        logger.error("No ORCA file has been provided")
        sys.exit(1)

    with open(filename, "r") as f:
        orca_file = [line.strip() for line in f]

    # Get the program version
    x = get_values(orca_file, "Program Version").result[0]
    x = x.replace("Program Version", "").strip()
    logger.info(f"Simulation run with ORCA version {x}")
    logger.info(mf.dashes)

    # Summary of the run
    run_summary(orca_file)

    if args["extract"] is not None:
        # Extract selected data
        extract_cmds = list_to_kwargs(args["extract"])

        for k, flags in extract_cmds.items():
            logger.info(mf.dashes)

            if k == "input":
                logger.info("Extracting Input Commands")
                input_commands = get_input_commands(orca_file)
                input_commands.write()

            elif k == "charges":
                logger.info("Extracting Atomic Charges")
                charges = get_charges(orca_file)
                charges.write()

            elif k == "coord":
                logger.info("Extracting Coordinates")
                coordinates = get_coordinates(orca_file)
                coordinates.write()

            elif k == "gap":
                logger.info("Extracting HOMO-LUMO Gap")
                orbital_energies = get_orbital_energies(orca_file)
                orbital_energies.write()

            elif k == "mo":
                logger.info("Extracting Molecular Orbitals")
                orbital_energies = get_orbital_energies(orca_file)
                orbital_energies.compute_mo(**flags)

            elif k == "nto":
                logger.info("Extracting Natural Transiotion Orbitals")
                tddft = get_tddft(orca_file)
                tddft.compute_nto(**flags)

            elif k == "animate":
                vibrational_modes = get_vibrational_modes(orca_file)
                vibrational_modes.animate(**flags)

            elif k == "+follow":
                vibrational_modes = get_vibrational_modes(orca_file)
                vibrational_modes.follow_negative_modes(**flags)
