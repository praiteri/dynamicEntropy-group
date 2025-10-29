#!/usr/bin/env python3
import io
import sys
import argparse
import logging
import pandas as pd


# Included from: from data_tools.ti_analysis.ti_analysis import ti_analysis
import sys
import logging
import numpy as np

from scipy.constants import hbar, electron_volt, atomic_mass, Boltzmann
from scipy.integrate import simpson

from chempy import Substance
from chempy.util.periodic import symbols


def get_dG_com(T, N, V, species, k):
    """
    Calculate center-of-mass free energy correction (δF_CM) in eV using mass fractions.

    Parameters:
    - T: temperature (K)
    - N: number of molecules
    - V: volume in Å³
    - species: list of species chempy format for one molecule
    - k: spring constant in eV/Å²

    Returns:
    - δF_CM in eV
    """
    # Create array with the massess of all atoms
    masses = []
    for atomic_number, count in species.composition.items():
        element_name = symbols[atomic_number - 1]
        element = Substance.from_formula(element_name)
        masses.extend([element.mass] * count)

    # Calculate total mass
    total_mass = sum(masses)

    # Calculate mass fractions (μᵢ = mᵢ/∑ᵢ₌₁ᴺmᵢ)
    mu_i = [m / total_mass for m in masses] * N

    # Calculate sum of squared mass fractions
    sum_mu_squared = sum(mu**2 for mu in mu_i)

    # Calculate according to equation
    k_SI = k * electron_volt * 1e20
    term = k_SI / (2 * np.pi * Boltzmann * T * sum_mu_squared)

    # Convert volume to SI units (from Å^3)
    V_m3 = V * 1e30
    prefactor = V_m3 / N
    delta_F_J = -Boltzmann * T * np.log(prefactor * term ** (3 / 2))

    return delta_F_J / electron_volt / N


def get_einstein_dG(species, k=1.0, T=300):
    """
    Calculate Einstein crystal free energy.
    """
    logger = logging.getLogger("main")
    # fmt = formatting().fmt

    # Assuming input k is in eV
    # Converting spring constant to SI units
    k_SI = k * electron_volt * 1e20  # J/m² (kg/s²)
    logger.debug(f"  Spring constant (eV): {k}")

    einstein_dG = 0.0
    for atomic_number, count in species.composition.items():
        element_name = symbols[atomic_number - 1]
        element = Substance.from_formula(element_name)

        logger.debug(f"  Number of {element_name} in formula unit: {count}")
        mass_kg = element.mass * atomic_mass
        omega = np.sqrt(k_SI / mass_kg)
        logger.debug(f"  Harmonic oscillator frequency (Hz): {omega:.3e}")

        dG = 3 * Boltzmann * T * np.log((hbar * omega) / (Boltzmann * T))
        logger.debug(f"  Free energy (eV): {dG / electron_volt:.5f}")
        dG *= count

        einstein_dG += dG

    return einstein_dG / electron_volt


def ti_analysis(df, **kwargs):
    # Initialise formatter
    logger = logging.getLogger("main")

    # Set parameters
    args = {
        "timestep": float(kwargs.get("ts", 0.001)),
        "temp": float(kwargs.get("temp", 300)),
        "kappa": float(kwargs.get("kappa", 1.0)),
        "nmols": int(kwargs.get("nmols", 1)),
        "chem": kwargs.get("chem"),
        "volume": float(kwargs.get("vol", 0)),
        "mode": kwargs.get("mode", "0"),
        "plot": kwargs.get("plot"),
    }

    logger.info(f"LAMMPS TI mode: {args['mode']}")
    logger.info(f"Timestep: {args['timestep']}")
    logger.info(f"Temperature (K): {args['temp']}")
    logger.info(f"Spring constant (eV): {args['kappa']}")
    logger.info(f"Volume (Å³): {args['volume']}")
    if args["nmols"] is None:
        logger.warning("Number of molecules not specified, assuming 1")
        args["nmols"] = 1
    else:
        args["nmols"] = int(args["nmols"])
    if args["chem"] is None:
        logger.warning(
            "Chemical formula not specified, no Einstein crystal free energy will be computed"
        )
    else:
        logger.info(f"Chemical formula: {args['chem']}")
        try:
            species = Substance.from_formula(args["chem"])
            string = ""
            for atomic_number, _ in species.composition.items():
                string += f"{symbols[atomic_number - 1]} "

            logger.info(f"Species: {string.strip()}")
        except Exception as e:
            logger.error(f"Error parsing chemical formula: {e}")
            sys.exit(1)

    # Chemical formula of the species
    species = None
    dg_einstein = 0.0
    if args["chem"] is not None:
        logger.info("Computing free energy of Einstein crystal")
        species = Substance.from_formula(args["chem"])

        # Compute harmonic oscillator free energies
        dg_einstein = get_einstein_dG(species, k=args["kappa"], T=args["temp"])
        logger.info(f"dG Einstein crystal : {dg_einstein:8.5f}")

    # Calculate COM correction using mass fractions
    dG_com_total = 0.0
    if args["volume"] > 0 and args["chem"] is not None:
        logger.info("Computing centre of mass correction term")
        dG_com_total = get_dG_com(
            args["temp"], args["nmols"], args["volume"], species, k=args["kappa"]
        )
        logger.info(f"Centre of mass correction : {dG_com_total:8.5f}")

    # Processing section
    step = df["Step"]
    time = df["Time"]
    erg_system = df["PotEng"]
    erg_springs = df["f_ti"]
    lambda1 = df["f_ti[1]"]
    dlambda1_ds = df["f_ti[2]"]

    if args["mode"] == "0":
        mask_fwd = dlambda1_ds > 0
        mask_bwd = dlambda1_ds < 0

        y = (erg_springs - erg_system) * dlambda1_ds / args["timestep"]
        fwd = simpson(y=y[mask_fwd], x=time[mask_fwd]) / args["nmols"]
        bwd = simpson(y=y[mask_bwd], x=time[mask_bwd]) / args["nmols"]

        logger.debug(f"  Forward TI (total) : {fwd:8.5f} (eV)")
        logger.debug(f"  Backward TI (total) : {bwd:8.5f} (eV)")
        dG_TI = 0.5 * (fwd - bwd)

    elif args["mode"] == "1":
        lambda1 = df["f_ti[1]"]
        dlambda1_ds = df["f_ti[2]"]
        lambda2 = df["f_ti[3]"]
        dlambda2_ds = df["f_ti[4]"]
        mask_fwd = (dlambda1_ds > 0) | (dlambda2_ds > 0)
        mask_bwd = (dlambda1_ds < 0) | (dlambda2_ds < 0)

        y_bulk = -erg_system * dlambda1_ds / args["timestep"]
        fwd_bulk = simpson(y=y_bulk[mask_fwd], x=time[mask_fwd]) / args["nmols"]
        bwd_bulk = simpson(y=y_bulk[mask_bwd], x=time[mask_bwd]) / args["nmols"]

        y_springs = erg_springs * dlambda2_ds / args["timestep"]
        fwd_springs = simpson(y=y_springs[mask_fwd], x=time[mask_fwd]) / args["nmols"]
        bwd_springs = simpson(y=y_springs[mask_bwd], x=time[mask_bwd]) / args["nmols"]

        logger.debug(f"  Forward TI (springs) : {fwd_springs:8.5f}")
        logger.debug(f"  Backward TI (springs) : {bwd_springs:8.5f}")
        logger.debug(f"  Forward TI (bulk) : {fwd_bulk:8.5f}")
        logger.debug(f"  Backward TI (bulk) : {bwd_bulk:8.5f}")
        logger.debug(f"  Forward TI (total) : {fwd_bulk + fwd_springs:8.5f}")
        logger.debug(f"  Backward TI (total) : {bwd_bulk + bwd_springs:8.5f}")

        rtmp = 0.5 * (fwd_springs - bwd_springs)
        logger.info(f"Average TI (springs) : {rtmp:8.5f}")
        rtmp = 0.5 * (fwd_bulk - bwd_bulk)
        logger.info(f"Average TI (bulk) : {rtmp:8.5f}")
        dG_TI = 0.5 * (fwd_springs + fwd_bulk - bwd_bulk - bwd_springs)

    elif args["mode"] == "2":
        lambda1 = df["f_ti[1]"]
        dlambda1_ds = -df["f_ti[2]"]
        lambda2 = df["f_ti[3]"]
        dlambda2_ds = df["f_ti[4]"]
        mask_fwd = (dlambda1_ds > 0) | (dlambda2_ds > 0)
        mask_bwd = (dlambda1_ds < 0) | (dlambda2_ds < 0)

        y_bulk = -erg_system * dlambda1_ds / args["timestep"]
        fwd_bulk = simpson(y=y_bulk[mask_fwd], x=time[mask_fwd]) / args["nmols"]
        bwd_bulk = simpson(y=y_bulk[mask_bwd], x=time[mask_bwd]) / args["nmols"]

        y_springs = erg_springs * dlambda2_ds / args["timestep"]
        fwd_springs = simpson(y=y_springs[mask_fwd], x=time[mask_fwd]) / args["nmols"]
        bwd_springs = simpson(y=y_springs[mask_bwd], x=time[mask_bwd]) / args["nmols"]

        logger.debug(f"  Forward TI (springs) : {fwd_springs:8.5f}")
        logger.debug(f"  Backward TI (springs) : {bwd_springs:8.5f}")
        logger.debug(f"  Forward TI (bulk) : {fwd_bulk:8.5f}")
        logger.debug(f"  Backward TI (bulk) : {bwd_bulk:8.5f}")
        logger.debug(f"  Forward TI (total) : {fwd_bulk + fwd_springs:8.5f}")
        logger.debug(f"  Backward TI (total) : {bwd_bulk + bwd_springs:8.5f}")

        rtmp = 0.5 * (fwd_springs - bwd_springs)
        logger.info(f"Average TI (springs) : {rtmp:8.5f}")
        rtmp = 0.5 * (fwd_bulk - bwd_bulk)
        logger.info(f"Average TI (bulk) : {rtmp:8.5f}")
        dG_TI = 0.5 * (fwd_springs + fwd_bulk - bwd_bulk - bwd_springs)

    else:
        logger.error(f"Unknown run mode {args['mode']}")
        sys.exit(1)

    logger.info(f"Average TI (total) : {dG_TI:8.5f} (eV)")
    total = dG_TI - dg_einstein + dG_com_total
    logger.info(f"Grand total : {total:8.5f}")
    logger.warning("All energies are in eV and per molecule")



def read_lammps_output_to_df(file_path):
    # Read the entire file
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Find all Step-Loop blocks
    data_blocks = []
    i = 0
    while i < len(lines):
        if lines[i].strip().startswith("Step"):
            header_index = i

            # Find the end of this block
            footer_index = -1
            for j in range(header_index + 1, len(lines)):
                if lines[j].strip().startswith("Loop") or lines[j].strip().startswith("WARNING"):
                    footer_index = j
                    break

            if footer_index == -1:
                # If no end marker found, use till the end of file
                data_lines = lines[header_index:]
                i = len(lines)  # Exit the loop after processing
            else:
                # Use lines between header and footer
                data_lines = lines[header_index:footer_index]
                i = footer_index + 1  # Move past this block

            # Convert the extracted portion into a DataFrame
            data_text = "".join(data_lines)
            try:
                df = pd.read_csv(io.StringIO(data_text), sep="\\s+")
                data_blocks.append(df)
            except Exception as e:
                try:
                    logger = logging.getLogger("main")
                    logger.info(f"Warning: Failed to parse a data block: {e}")
                except NameError:
                    print(f"Warning: Failed to parse a data block: {e}")
        else:
            i += 1

    if not data_blocks:
        raise ValueError("Could not find any lines starting with 'Step'")

    if len(data_blocks) == 1:
        data_blocks = data_blocks[0]

    return data_blocks


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--logfile", "--log", type=str, default=None, help="Log file")
    parser.add_argument("--quiet", "--q", action="store_true", help="Quiet mode")
    parser.add_argument("--debug", "--vv", action="store_true", help="Debug mode")

    parser.add_argument("--temp", "--t", type=float, help="Temperature")
    parser.add_argument("--kappa", "--k", type=float, help="Spring constant eV/Å^2")

    parser.add_argument("--input", "--i", type=str, help="Lammps output files with the TI data")
    parser.add_argument("--mode", "--m", type=str, help="TI simulation mode")
    parser.add_argument("--timestep", "--ts", type=float, help="Timestep")

    parser.add_argument("--nmols", "--n", type=int, help="Number of molecules")
    parser.add_argument("--chem", "--s", type=str, help="Chemical formula")
    parser.add_argument("--volume", "--vol", type=float, help="Volume of the simulation")

    args = parser.parse_args()

    logger = logging.getLogger("main")
    if args.quiet:
        logger.setLevel(logging.ERROR)
    elif args.debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    formatter = logging.Formatter("%(levelname)s - %(message)s")
    if args.logfile:
        handler = logging.FileHandler(args.logfile, mode='w')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    else:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    if args.input is None:
        logger.error("No input file provided. Use --input to specify the LAMMPS output file.")
        sys.exit(1)

    logger.info(f"Reading data from {args.input}")
    df = read_lammps_output_to_df(args.input)
    logger.info("LAMMPS output file read successfully.")
    logger.info(f"\n{df.head()}\n...")
    # Here you would call the ti_analysis function with the parsed arguments
    args_dict = {k: v for k, v in vars(args).items() if v is not None}

    ti_analysis(df, **args_dict)


if __name__ == "__main__":
    main()
