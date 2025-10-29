import argparse

# from ..orca_tools.orca_config import OrcaConfig


def setup_argument_parser():
    """Set up and return the argument parser for Orca calculations."""
    parser = argparse.ArgumentParser(description="Run Orca calculations")

    # Create argument groups
    mode_group = parser.add_argument_group(
        "Mode Options", "Debug and verbosity options"
    )
    # orca_config_group = parser.add_argument_group(
    #     "Orca Configuration", "Configuration parameters for Orca calculations"
    # )
    # run_options_group = parser.add_argument_group(
    #     "Run Options", "Workflow and execution options"
    # )
    post_processing_group = parser.add_argument_group(
        "Post-processing", "Options for analyzing output files"
    )

    # # Get argument configurations from OrcaConfig and add them to the orca_config_group
    # arg_info = OrcaConfig.get_arg_info()
    # for param_name, config in arg_info.items():
    #     orca_config_group.add_argument(f"--{param_name}", **config)

    # Add arguments to the mode_group
    mode_group.add_argument("--debug", action="store_true", help="Debug mode")
    mode_group.add_argument(
        "--verbose", "--v", action="store_true", help="Verbose mode"
    )

    # Add arguments to the post_processing_group
    post_processing_group.add_argument(
        "--input", "--i", type=str, help="Analyse ORCA output file"
    )
    post_processing_group.add_argument(
        "--extract", type=str, help="Extract data from ORCA output file", nargs="*"
    )

    # Add arguments to the run_options_group
    # run_options_group.add_argument(
    #     "--run", action="store_true", help="Create ORCA input"
    # )
    # run_options_group.add_argument(
    #     "--dry", action="store_true", help="Dry run (do not run ORCA)"
    # )
    # run_options_group.add_argument(
    #     "--wf", type=str, default=None, help="Calculation workflow"
    # )

    args = parser.parse_args()

    return args
