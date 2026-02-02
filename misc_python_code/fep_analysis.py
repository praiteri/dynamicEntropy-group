#!/usr/bin/env python3
"""
Free Energy Perturbation (FEP) Analysis Tool

This script performs BAR (Bennett Acceptance Ratio) and MBAR (Multistate Bennett
Acceptance Ratio) analysis on FEP simulation data across multiple lambda windows.
Supports analysis of different FEP types including coulombic and van der Waals interactions.

Custom FEP Folders:
  - By default, the script analyzes standard FEP stages (FEP_vdw, FEP_coul, etc.)
  - Custom folders can be added via command line using --custom-folder or --custom-folder-sign
  - Custom folders default to sign=1 unless specified otherwise
  - Use --only-custom to analyze only custom folders (ignoring defaults)
  
Programmatic Usage:
  - add_fep_stage(folder_name, sign=1): Add a custom FEP folder
  - set_fep_stages(stages_dict): Replace entire FEP stages configuration
  - reset_fep_stages(): Reset to default configuration
  - get_fep_stages(): Get current FEP stages configuration
"""

import os
import sys
import glob
import argparse
import logging
from typing import Dict, List, Tuple, Optional, Callable

import numpy as np
import pandas as pd
import pymbar as mb
from colorama import Fore, Style


# ============================================================================
# Configuration Constants
# ============================================================================

INDENTATION = 2
LEVEL_NAME_WIDTH = 8

# Default FEP stages with their sign conventions
FEP_STAGES = {
    "FEP_vdw": 1,
    "FEP_coul": 1,
    "FEP_coul_vacuum": -1,
    "FEP_file": 1,
    "FEP_ff": 1,
}


def get_fep_stages() -> Dict[str, int]:
    """
    Get the current FEP stages configuration.
    
    Returns:
        Dictionary mapping folder names to their sign conventions
    """
    return FEP_STAGES.copy()


def add_fep_stage(folder_name: str, sign: int = 1) -> None:
    """
    Add a custom FEP stage to the analysis.
    
    Args:
        folder_name: Name of the FEP folder to add
        sign: Sign convention for the stage (default: 1)
    """
    FEP_STAGES[folder_name] = sign
    logger = logging.getLogger("fep_analysis")
    logger.info(f"Added custom FEP stage: {folder_name} (sign={sign})")


def set_fep_stages(stages: Dict[str, int]) -> None:
    """
    Replace the entire FEP stages configuration.
    
    Args:
        stages: Dictionary mapping folder names to their sign conventions
    """
    global FEP_STAGES
    FEP_STAGES = stages.copy()
    logger = logging.getLogger("fep_analysis")
    logger.info(f"Updated FEP stages configuration: {list(stages.keys())}")


def reset_fep_stages() -> None:
    """Reset FEP stages to default configuration."""
    global FEP_STAGES
    FEP_STAGES = {
        "FEP_vdw": 1,
        "FEP_coul": 1,
        "FEP_coul_vacuum": -1,
        "FEP_file": 1,
        "FEP_ff": 1,
    }
    logger = logging.getLogger("fep_analysis")
    logger.info("Reset FEP stages to default configuration")


UNIT_SYSTEMS = {
    "omm": {
        "energy": "kJ/mol",
        "distance": "nm",
        "c0": 0.6022,  # 1/(nm**3)
        "kB": 0.00831446261815324,  # kJ/(K mol)
        "pi4e0": 138.93545302827604,  # nm kJ/(mol e**2)
    },
    "lmp": {
        "energy": "eV",
        "distance": "A",
        "c0": 0.0006022,  # 1/(A**3)
        "kB": 0.000086173,  # eV / K
        "pi4e0": 14.399645,  # A eV/(e**2)
    },
}

DIELECTRIC_MODELS = {
    "exp": lambda T: 249.355 - 0.7877 * T + 0.0007192 * T**2,
    "spcfw": lambda T: 181.464 - 0.338253 * T,
    "amoeba": lambda T: 412.207 - 1.5853 * T + 0.0017023 * T**2,
}


# ============================================================================
# Logging Configuration
# ============================================================================


class ColoredFormatter(logging.Formatter):
    """Custom formatter that adds colors to log messages based on level."""

    COLORS = {
        "DEBUG": Fore.CYAN + Style.BRIGHT,
        "INFO": Fore.GREEN + Style.BRIGHT,
        "WARNING": Fore.YELLOW + Style.BRIGHT,
        "ERROR": Fore.RED + Style.BRIGHT,
        "CRITICAL": Fore.MAGENTA + Style.BRIGHT,
    }

    def format(self, record: logging.LogRecord) -> str:
        """Format log record with color coding."""
        orig_levelname = record.levelname
        orig_msg = record.msg

        color_code = self.COLORS.get(orig_levelname, "")
        padding = " " * (LEVEL_NAME_WIDTH - len(orig_levelname))

        record.levelname = f"{color_code}{orig_levelname}{padding}{Style.RESET_ALL}"
        record.msg = f"{color_code}{orig_msg}{Style.RESET_ALL}"

        result = super().format(record)

        record.levelname = orig_levelname
        record.msg = orig_msg
        return result


def setup_logger(debug: bool = False) -> logging.Logger:
    """
    Configure and return application logger.

    Args:
        debug: Enable debug-level logging if True

    Returns:
        Configured logger instance
    """
    logger = logging.getLogger("fep_analysis")
    logger.setLevel(logging.DEBUG if debug else logging.INFO)
    logger.handlers.clear()

    handler = logging.StreamHandler(stream=sys.stdout)
    formatter = ColoredFormatter(fmt="%(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    return logger


# ============================================================================
# Logging Utilities
# ============================================================================


def log_message(
    data: Optional[any] = None,
    level: str = "INFO",
    indent: int = 0,
    title: Optional[str] = None,
    marker: str = ":",
    align_width: Optional[int] = None,
    separator: bool = False,
    ncols: int = 6,
) -> None:
    """
    Pretty print data with optional title and formatting.

    Args:
        data: Data to log (dict, list, or scalar)
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        indent: Indentation level
        title: Optional title to display
        marker: Separator between key and value
        align_width: Column width for alignment
        separator: Print separator line if True
        ncols: Number of columns for list display
    """
    logger = logging.getLogger("fep_analysis")
    log_func = getattr(logger, level.lower(), logger.info)

    if separator:
        log_func("=" * 60)

    if title is not None:
        indent_str = " " * (indent * INDENTATION)
        log_func(f"{indent_str}{title}")

    if data is None:
        return

    if isinstance(data, dict):
        _log_dict(
            data, log_func, indent + (1 if title else 0), marker, align_width, ncols
        )
    elif isinstance(data, list):
        _log_list(data, log_func, indent + (1 if title else 0), ncols)
    else:
        indent_str = " " * (indent * INDENTATION)
        log_func(f"{indent_str}{data}")


def _log_dict(
    data: Dict,
    log_func: Callable,
    indent: int,
    marker: str,
    align_width: Optional[int],
    ncols: int,
) -> None:
    """Log dictionary with proper formatting."""
    indent_str = " " * (indent * INDENTATION)
    align_col = (align_width or 48) - len(indent_str)

    for key, value in data.items():
        if isinstance(value, dict):
            log_func(f"{indent_str}{key}{marker}")
            _log_dict(value, log_func, indent + 1, marker, align_width, ncols)
        elif isinstance(value, list):
            log_func(f"{indent_str}{key}{marker}")
            _log_list(value, log_func, indent + 1, ncols)
        elif isinstance(value, bool):
            log_func(f"{indent_str}{key:<{align_col}} {marker} {str(value).lower()}")
        elif value is None:
            log_func(f"{indent_str}{key:<{align_col}} {marker} null")
        else:
            log_func(f"{indent_str}{key:<{align_col}} {marker} {value}")


def _log_list(data: List, log_func: Callable, indent: int, ncols: int) -> None:
    """Log list in columnar format."""
    indent_str = " " * (indent * INDENTATION)
    for i in range(0, len(data), ncols):
        chunk = data[i : i + ncols]
        string = " ".join(f"{x:8}" for x in chunk)
        log_func(f"{indent_str}{string}")


# ============================================================================
# Data I/O
# ============================================================================


def read_mbar_data(filepath: str) -> Tuple[np.ndarray, List[str], pd.DataFrame]:
    """
    Read FEP output file and extract energy data.

    Args:
        filepath: Path to FEP output file

    Returns:
        Tuple of (energies array, lambda column names, full dataframe)

    Raises:
        ValueError: If no energy columns found in file
    """
    df = pd.read_csv(filepath, sep=r"\s+", comment="#")
    lambda_cols = [col for col in df.columns if col.startswith("E")]

    if not lambda_cols:
        raise ValueError(f"No energy columns found in {filepath}")

    energies = df[lambda_cols].values
    return energies, lambda_cols, df


def read_bar_data(filepath: str) -> pd.DataFrame:
    """
    Read and parse BAR FEP output file.

    Args:
        filepath: Path to FEP output file

    Returns:
        DataFrame with time, forward, and backward energy differences
    """
    df = pd.read_csv(filepath, sep=r"\s+")

    # Handle column shift issue from merged data
    if df.columns[0] == "#####":
        tmp = df.columns[1:]
        df = df.drop(df.columns[-1], axis=1)
        df.columns = tmp

    # Try different column name conventions
    try:
        return pd.DataFrame(
            {
                "time": df["Time"],
                "du_fwd": df["E1"] - df["E0"],
                "du_bwd": df["E2"] - df["E0"],
            }
        )
    except KeyError:
        return pd.DataFrame(
            {
                "time": df["Time"],
                "du_fwd": df["Ef"] - df["E0"],
                "du_bwd": df["Eb"] - df["E0"],
            }
        )


def read_fep_data(filepath: str) -> pd.DataFrame:
    """
    Read and parse BAR FEP output file.

    Args:
        filepath: Path to FEP output file

    Returns:
        DataFrame with time, forward, and backward energy differences
    """
    df = pd.read_csv(filepath, sep=r"\s+", comment="#")
    lambda_cols = [col for col in df.columns if col.startswith("E")]
    result = pd.DataFrame(
        {
            "time": df["Time"],
        }
    )
    for c in lambda_cols[1:]:
        new_col = f"({c}-{lambda_cols[0]})"
        result[new_col] = df[c] - df[lambda_cols[0]]
    return result


# ============================================================================
# Analysis Functions
# ============================================================================


def mbar_analysis(params: Dict) -> Dict:
    """
    Perform MBAR (Multistate Bennett Acceptance Ratio) analysis.

    Args:
        params: Dictionary containing analysis parameters including:
            - units: Unit system specification
            - temp: Temperature
            - equil: Equilibration time
            - skip: Frame skip interval

    Returns:
        Dictionary with stage names, free energies, and uncertainties
    """
    kT = params["units"]["kB"] * params["temp"]
    beta = 1.0 / kT

    results = {"stage": [], "ΔG": [], "err": []}

    log_message(separator=True)
    for folder, sign in FEP_STAGES.items():
        if not os.path.isdir(folder):
            continue

        # Find MBAR output files
        fep_files = sorted(glob.glob(os.path.join(folder, "mbar.?.out")))
        fep_files.extend(sorted(glob.glob(os.path.join(folder, "mbar.??.out"))))
        n_states = len(fep_files)

        if n_states == 0:
            continue

        log_message(f"Found {n_states} lambda states in {folder}", level="debug")

        # Read and process energy data
        u_kn_list = []
        N_k = np.zeros(n_states, dtype=int)

        for lambda_num, fep_file in enumerate(fep_files):
            if not os.path.exists(fep_file):
                raise ValueError(f"File not found: {fep_file}")

            energies, _, df = read_mbar_data(fep_file)

            if energies.shape[1] != n_states:
                raise ValueError(
                    f"Expected {n_states} energy columns, got {energies.shape[1]}"
                )

            # Apply equilibration and skip filters
            mask = df["Time"] > params["equil"]
            energies_processed = energies[mask][:: params["skip"]]
            energies_processed = energies_processed.T * beta

            u_kn_list.append(energies_processed)
            N_k[lambda_num] = energies_processed.shape[1]

        # Perform MBAR analysis
        u_kn = np.concatenate(u_kn_list, axis=1)
        mbar = mb.MBAR(u_kn, N_k, verbose=True, initialize="BAR")
        mbar_results = mbar.compute_free_energy_differences()

        delta_f = mbar_results["Delta_f"]
        d_delta_f = mbar_results["dDelta_f"]

        # Convert to energy units
        delta_f_energy = delta_f * kT * sign
        d_delta_f_energy = d_delta_f * kT

        # Store results
        results["stage"].append(folder)
        results["ΔG"].append(delta_f_energy[0, -1])
        results["err"].append(d_delta_f_energy[0, -1])

        # Display results
        _display_mbar_results(
            folder, n_states, delta_f_energy, d_delta_f_energy, params
        )

    return results


def _display_mbar_results(
    folder: str,
    n_states: int,
    delta_f: np.ndarray,
    d_delta_f: np.ndarray,
    params: Dict,
) -> None:
    """Display MBAR analysis results."""
    # log_message(separator=True)
    # # log_message(f"MBAR RESULTS for {folder}")
    # # log_message(separator=True)
    log_message(
        f"Free energy differences ({params['units']['energy']}) relative to state 0:",
        level="debug",
    )
    log_message(f"{'State':<5} {'ΔF':>15}  {'ΔΔF':>10} {'dΔF':>15}", level="debug")
    log_message("-" * 50, level="debug")
    for i in range(n_states):
        log_message(
            f"{i:>5} {delta_f[0, i]:>15.4f} [{delta_f[0, i]-delta_f[0, max(0,i-1)]:>10.4f}] {d_delta_f[0, i]:>15.4f}",
            level="debug",
        )

    log_message(
        {
            f"Total for {folder:15}": f"{delta_f[0, -1]:8.3f} +/- {d_delta_f[0, -1]:8.3f}"
        },
    )


def bar_analysis(params: Dict) -> Dict:
    """
    Perform BAR (Bennett Acceptance Ratio) analysis.

    Args:
        params: Dictionary containing analysis parameters

    Returns:
        Dictionary with stage names, free energies, and uncertainties
    """
    kT = params["units"]["kB"] * params["temp"]
    beta = 1.0 / kT

    results = {
        "stage": [],
        "ΔG": [],
        "ΔG_f": [],
        "ΔG_b": [],
        "err": [],
        "err_f": [],
        "err_b": [],
    }

    log_message(separator=True)
    for folder, sign in FEP_STAGES.items():
        if not os.path.isdir(folder):
            continue

        # Find FEP output files
        fep_files = sorted(glob.glob(os.path.join(folder, "bar.?.out")))
        fep_files.extend(sorted(glob.glob(os.path.join(folder, "bar.??.out"))))

        if not fep_files:
            continue

        # Read all stage data
        stage_data = []
        for i, filepath in enumerate(fep_files):
            df = read_bar_data(filepath)
            df["du_fwd"] *= beta
            df["du_bwd"] *= beta
            stage_data.append(df)
            results["stage"].append(filepath)

        # Analyze each transition
        total_dg = 0.0
        total_err_sq = 0.0

        for i in range(len(stage_data) - 1):
            # Extract forward and backward work values
            df_current = stage_data[i]
            df_next = stage_data[i + 1]

            w_F = df_current[df_current["time"] > params["equil"]]["du_fwd"][
                :: params["skip"]
            ]
            w_B = df_next[df_next["time"] > params["equil"]]["du_bwd"][
                :: params["skip"]
            ]

            # Compute estimates
            fwd_est = mb.other_estimators.exp(w_F, is_timeseries=True)
            bar_est = mb.other_estimators.bar(w_F, w_B, uncertainty_method="BAR")
            bwd_est = mb.other_estimators.exp(w_B, is_timeseries=True)

            # Store results
            results["ΔG"].append(kT * bar_est["Delta_f"] * sign)
            results["ΔG_f"].append(kT * fwd_est["Delta_f"] * sign)
            results["ΔG_b"].append(kT * bwd_est["Delta_f"] * sign)
            results["err"].append(kT * bar_est["dDelta_f"])
            results["err_f"].append(kT * fwd_est["dDelta_f"])
            results["err_b"].append(kT * bwd_est["dDelta_f"])

            # Display individual transition
            filepath = fep_files[i]
            log_message(
                {
                    f"File {filepath:20}": f"{results['ΔG'][-1]:8.3f} +/- {results['err'][-1]:8.3f}"
                },
                level="debug",
            )

            total_dg += results["ΔG"][-1]
            total_err_sq += results["err"][-1] ** 2

        # Display folder summary
        log_message(
            {
                f"Total for {folder:15}": f"{total_dg:8.3f} +/- {np.sqrt(total_err_sq):8.3f}"
            },
        )

    return results


def fep_analysis(params: Dict, arg) -> Dict:
    """
    Perform BAR (Bennett Acceptance Ratio) analysis.

    Args:
        params: Dictionary containing analysis parameters

    Returns:
        Dictionary with stage names, free energies, and uncertainties
    """
    kT = params["units"]["kB"] * params["temp"]
    beta = 1.0 / kT

    results = {
        "stage": [],
        "ΔG": [],
        "err": [],
    }

    log_message(separator=True)
    for folder, sign in FEP_STAGES.items():
        if not os.path.isdir(folder):
            continue

        # Find FEP output files
        fep_files = sorted(glob.glob(os.path.join(folder, "fep.?.out")))
        fep_files.extend(sorted(glob.glob(os.path.join(folder, "fep.??.out"))))

        if not fep_files:
            continue

        # Read all stage data
        stage_data = []
        for i, filepath in enumerate(fep_files):
            df = read_fep_data(filepath)
            for i in range(1, len(df.columns)):
                df.iloc[:, i] *= beta
            stage_data.append(df)
            results["stage"].append(filepath)

        # Analyze each transition
        total_dg = 0.0
        total_err_sq = 0.0

        index = None
        if isinstance(arg, bool) and arg:
            index = 1
        elif isinstance(arg, int):
            index = arg
        elif isinstance(arg, str):
            index = list(df.columns).index(arg)

        if index is None or index > (len(df.columns)):
            log_message(f"Cannot compute FEP for column {index}", level="error")
            sys.exit(1)

        c = df.columns[index]
        for i in range(len(stage_data)):
            # Extract forward and backward work values
            df_current = stage_data[i]

            w_F = df_current[df_current["time"] > params["equil"]][c][:: params["skip"]]
            if np.mean(w_F) == 0 and np.std(w_F) == 0:
                results["ΔG"].append(0)
                results["err"].append(0)
                continue

            # Compute estimates
            fwd_est = mb.other_estimators.exp(w_F, is_timeseries=True)

            # Store results
            results["ΔG"].append(kT * fwd_est["Delta_f"] * sign)
            results["err"].append(kT * fwd_est["dDelta_f"])

            # Display individual transition
            filepath = fep_files[i]
            log_message(
                {
                    f"File {filepath:20}": f"{results['ΔG'][-1]:8.3f} +/- {results['err'][-1]:8.3f}",
                },
                level="debug",
            )

            total_dg += results["ΔG"][-1]
            total_err_sq += results["err"][-1] ** 2

        # Display folder summary
        log_message(
            {
                f"Total for {folder:15}": f"{total_dg:8.3f} +/- {np.sqrt(total_err_sq):8.3f}"
            }
        )

    return results


# ============================================================================
# Parameter Configuration
# ============================================================================


def calculate_dielectric_constant(temperature: float, water_model: str) -> float:
    """
    Calculate dielectric constant for water model at given temperature.

    Args:
        temperature: Temperature in Kelvin
        water_model: Water model name (exp, spcfw, amoeba) or numerical value

    Returns:
        Dielectric constant value

    Raises:
        ValueError: If water model is not recognized
    """
    # Try to evaluate as numeric expression
    try:
        return float(eval(str(water_model)))
    except (NameError, SyntaxError, TypeError):
        pass

    # Look up model
    if water_model not in DIELECTRIC_MODELS:
        available = ", ".join(DIELECTRIC_MODELS.keys())
        raise ValueError(
            f"Unknown water model '{water_model}'. Available models: {available}"
        )

    return DIELECTRIC_MODELS[water_model](temperature)


def configure_parameters(args: argparse.Namespace) -> Dict:
    """
    Configure analysis parameters from command line arguments.

    Args:
        args: Parsed command line arguments

    Returns:
        Dictionary of analysis parameters
    """
    # Default parameters
    params = {
        "units": "omm",
        "temp": 300,
        "charge": 0.0,
        "water": "amoeba",
        "cell": None,
        "skip": 1,
        "equil": 0,
        "debug": False,
        "zeta": 2.837297,
    }

    # Override with command line arguments
    for key, value in vars(args).items():
        if value is not None and key in params:
            params[key] = value

    # Configure unit system
    if params["units"] not in UNIT_SYSTEMS:
        available = ", ".join(UNIT_SYSTEMS.keys())
        raise ValueError(
            f"Unknown unit system '{params['units']}'. Available: {available}"
        )

    params["units"] = UNIT_SYSTEMS[params["units"]]

    # Calculate dielectric constant
    params["eps"] = calculate_dielectric_constant(params["temp"], params["water"])

    log_message(params, level="INFO", title="Analysis Parameters")
    return params


# ============================================================================
# Command Line Interface
# ============================================================================


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.

    Returns:
        Parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Free Energy Perturbation (FEP) Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic BAR analysis
  %(prog)s --bar --temp 300 --skip 10 --equil 100
  
  # MBAR analysis with custom unit system
  %(prog)s --mbar --units lmp --water spcfw
  
  # Add custom FEP folders (sign=1 by default)
  %(prog)s --bar --custom-folder FEP_custom1 --custom-folder FEP_custom2
  
  # Add custom folder with specific sign
  %(prog)s --bar --custom-folder-sign FEP_custom -1
  
  # Analyze only custom folders (ignore defaults)
  %(prog)s --bar --only-custom --custom-folder FEP_my_analysis
        """,
    )

    def int_or_str(value):
        try:
            return int(value)
        except ValueError:
            return value

    # System parameters
    parser.add_argument("--charge", "-q", type=float, help="System charge")
    parser.add_argument(
        "--temperature",
        "-t",
        dest="temp",
        type=float,
        help="Temperature in Kelvin (default: 300)",
    )
    parser.add_argument(
        "--water",
        "-w",
        type=str,
        help="Water model (exp, spcfw, amoeba) or custom dielectric constant",
    )
    parser.add_argument("--cell", "-l", type=float, help="Cell size")
    parser.add_argument(
        "--units",
        "-u",
        choices=list(UNIT_SYSTEMS.keys()),
        help="Unit system (default: omm)",
    )

    # Analysis parameters
    parser.add_argument(
        "--skip", "-s", type=int, help="Frame skip interval (default: 1)"
    )
    parser.add_argument(
        "--equil", "-e", type=float, help="Equilibration time to discard (default: 0)"
    )
    
    # FEP stage configuration
    parser.add_argument(
        "--custom-folder",
        "-c",
        action="append",
        metavar="FOLDER",
        help="Add a custom FEP folder to analyze (sign=1). Can be specified multiple times."
    )
    parser.add_argument(
        "--custom-folder-sign",
        action="append",
        nargs=2,
        metavar=("FOLDER", "SIGN"),
        help="Add a custom FEP folder with specific sign. Format: FOLDER SIGN (e.g., 'FEP_custom -1')"
    )
    parser.add_argument(
        "--only-custom",
        action="store_true",
        help="Only analyze custom folders (ignore default FEP_STAGES)"
    )

    # Analysis method
    analysis_group = parser.add_mutually_exclusive_group(required=True)
    analysis_group.add_argument(
        "--fep",
        dest="run_fep",
        nargs="?",
        const=True,
        default=False,
        type=int_or_str,
        help="Run FEP analysis (optionally specify an integer or string)",
    )
    analysis_group.add_argument(
        "--bar",
        dest="run_bar",
        action="store_true",
        help="Run BAR (Bennett Acceptance Ratio) analysis",
    )
    analysis_group.add_argument(
        "--mbar",
        dest="run_mbar",
        action="store_true",
        help="Run MBAR (Multistate Bennett Acceptance Ratio) analysis",
    )

    # Debugging
    parser.add_argument(
        "--debug", "-v", action="store_true", help="Enable debug output"
    )

    return parser.parse_args()


# ============================================================================
# Main Execution
# ============================================================================


def main() -> None:
    """Main execution function."""
    args = parse_arguments()
    logger = setup_logger(args.debug)

    logger.info("FEP Analysis Tool")
    log_message(separator=True)

    if args.debug:
        log_message(vars(args), level="DEBUG", title="Command Line Arguments", indent=1)
        log_message(separator=True)

    # Configure FEP stages based on command line arguments
    if args.only_custom:
        # Reset to empty if only analyzing custom folders
        set_fep_stages({})
    
    # Add custom folders with default sign (1)
    if args.custom_folder:
        for folder in args.custom_folder:
            add_fep_stage(folder, sign=1)
    
    # Add custom folders with specific signs
    if args.custom_folder_sign:
        for folder, sign_str in args.custom_folder_sign:
            try:
                sign = int(sign_str)
                if sign not in [-1, 1]:
                    logger.error(f"Sign must be -1 or 1, got {sign}")
                    return
                add_fep_stage(folder, sign=sign)
            except ValueError:
                logger.error(f"Invalid sign value: {sign_str}")
                return

    # Configure parameters
    params = configure_parameters(args)

    # Run analysis
    if args.run_bar:
        results = bar_analysis(params)
    elif args.run_mbar:
        results = mbar_analysis(params)
    elif args.run_fep:
        results = fep_analysis(params, args.run_fep)
    else:
        logger.error("No analysis method specified")
        return

    # Calculate and display total free energy
    total_dg = np.sum(results["ΔG"])
    total_err = np.sqrt(np.sum(np.array(results["err"]) ** 2))

    log_message(separator=True)
    log_message(
        f"Total Free Energy ({params['units']['energy']}): "
        f"{total_dg:10.3f} ± {total_err:10.3f}"
    )
    log_message(separator=True)


if __name__ == "__main__":
    main()
