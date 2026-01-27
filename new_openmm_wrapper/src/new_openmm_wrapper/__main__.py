### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import sys
import numpy as np
import yaml
import new_openmm_wrapper as my

import time


def format_time(seconds):
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = int(seconds % 60)
    milliseconds = int((seconds % 1) * 1000)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}.{milliseconds:03d}"


def main():
    """Main execution function."""

    # Start timer
    start = time.time()

    # Parse command line arguments
    args = my.parse_arguments()

    # Initialise logger
    if args.debug:
        logger = my.create_logger(level="DEBUG", filename=args.logfile)
    else:
        logger = my.create_logger(filename=args.logfile)

    # my.test_pretty_log()

    my.pretty_log(title="runOpenMM_lite", sep=True)
    my.pretty_log(vars(args), logger="DEBUG", title="Command line arguments:", sep=True)

    # Read input
    try:
        # Define the constructor function for the yaml loader
        def arange_constructor(loader, node):
            params = loader.construct_sequence(node)
            array = np.arange(params[0], params[1] + params[2] / 2, params[2]).tolist()
            return [float(f"{x:.3f}") for x in array]

        # Register it BEFORE loading the YAML
        yaml.add_constructor("!arange", arange_constructor, Loader=yaml.SafeLoader)

        with open(args.input, "r") as f:
            input_commands = yaml.safe_load(f)
    except Exception as e:
        print(e)
        logger.error(e)
        sys.exit(1)

    my.pretty_log(
        input_commands, logger="DEBUG", title=f"Input file ({args.input}):", sep=True
    )
    input_commands["debug"] = args.debug

    run_type = next(
        (
            rtype
            for rtype in ["energy", "fep", "md"]
            if rtype in input_commands and input_commands[rtype]
        ),
        None,
    )

    if run_type is None:
        setup = my.simulationSetup(input_commands)
        my.createSystem(setup)

    elif run_type == "energy":
        my.single_point_energy(input_commands)

    elif run_type == "md":
        my.molecular_dynamics(input_commands)

    elif run_type == "fep":
        my.compute_fep(input_commands)

    # Stop timer
    end = time.time()
    elapsed = end - start
    formatted = format_time(elapsed)

    my.pretty_log(None, title="Normal termination", logger="INFO", sep=True)
    my.pretty_log({"Execution time": formatted}, align_width=0)
    my.pretty_log(sep=True)
    sys.exit(0)


if __name__ == "__main__":
    main()
