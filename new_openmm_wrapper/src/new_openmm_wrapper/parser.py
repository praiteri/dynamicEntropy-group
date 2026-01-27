### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import argparse

def parse_arguments():
    """
    Parse command line arguments.

    Returns:
        argparse.Namespace: Parsed command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Generate FEP configuration files for multiple lambda windows",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        dest="input",
        type=str,
        default="input.yaml",
        help="Input YAML configuration file (default: mdp.yaml)",
    )

    parser.add_argument(
        "-v", "--debug", dest="debug", action="store_true", help="Enable debug logging"
    )

    parser.add_argument(
        "-l", "--log", dest="logfile", type=str, default=None, help="Redirect screen output to file"
    )

    args = parser.parse_args()



    return args

