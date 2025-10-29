import os
import sys
from pathlib import Path
import argparse


# ============================================================================
# Command-Line Interface
# ============================================================================
class FileExtensionAction(argparse.Action):
    """
    Custom argparse action to handle file extension extraction.

    This action can either:
    1. Extract the extension from the provided filename (for --i flag)
    2. Use the extension specified in the flag itself (for --iXXX flags)
    """

    def __init__(self, option_strings, dest, **kwargs):
        # Determine if we should use the file extension from the filename
        # or from the flag itself based on the option string
        if "--i" in option_strings[0] and len(option_strings[0]) != 3:
            self.use_file_extension = False  # Use extension from flag (--iXXX)
        else:
            self.use_file_extension = True  # Use extension from filename (--i)
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        if self.use_file_extension:
            # For --i command: extract extension from the filename
            _, ext = os.path.splitext(values[0])
            if not ext:
                parser.error(f"No extension found in filename: {values[0]}")
            ext = ext.lstrip(".")
        else:
            # For --iXXX commands: extract extension from the flag itself
            if "--i" in option_string:
                ext = option_string[3:]  # Remove '--i' prefix
            elif "-i" in option_string:
                ext = option_string[2:]  # Remove '-i' prefix
            else:
                parser.error(f"Invalid option string: {option_string}")
                sys.exit(1)

            if not ext:
                parser.error("No extension specified in flag")

        # Store both filename and extracted extension
        result = {
            "filename": values[0],
            "type": ext,
        }

        setattr(namespace, self.dest, result)


def basic_parser(description):
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("--quiet", action="store_true", help="Quiet mode")
    parser.add_argument("--debug", "--vv", action="store_true", help="Debug mode")
    parser.add_argument("--verbose", "--v", action="store_true", help="Verbose mode")
    parser.add_argument("--result", action="store_true", help="Result mode")

    for ext in ["", "txt", "dat", "csv", "lmp"]:  # Add more extensions as needed
        if len(ext) == 0:
            msg = "Input file (type is auto-detected from extension)"
        else:
            msg = f"Input file (forced {ext} extension)"
        parser.add_argument(
            f"--i{ext}",
            f"-i{ext}",
            action=FileExtensionAction,
            type=str,
            nargs="+",
            dest="input",
            help=msg,
        )

    parser.add_argument("--output", "--o", type=str, default=None, help="Quiet mode")

    # Column selection
    parser.add_argument(
        "--columns",
        "--columns",
        "--c",
        type=str,
        default=None,
        nargs="*",
        help="Specify columns to read (e.g., '1,2,5' or 'colA,colB' or '1 2' or 'Time Lx,Lz' etc.)",
    )
    parser.add_argument(
        "--select",
        "--sel",
        type=str,
        nargs="+",
        action="extend",
        help="Range selection for a column",
    )
    parser.add_argument(
        "--plottype",
        "--pt",
        type=str,
        default=None,
        choices=["plotly", "pl", "matplotlib", "pyplot", "pp"],
        help="Range selection for a column",
    )

    return parser


def basic_parser_2(description):
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("--quiet", action="store_true", help="Quiet mode")
    parser.add_argument("--debug", "--vv", action="store_true", help="Debug mode")
    parser.add_argument("--verbose", "--v", action="store_true", help="Verbose mode")
    parser.add_argument("--result", action="store_true", help="Result mode")

    for ext in ["", "txt", "dat", "csv", "lmp"]:  # Add more extensions as needed
        if len(ext) == 0:
            msg = "Input file (type is auto-detected from extension)"
        else:
            msg = f"Input file (forced {ext} extension)"
        parser.add_argument(
            f"--i{ext}",
            f"-i{ext}",
            # action=FileExtensionAction,
            type=str,
            nargs="+",
            dest=f"input_{ext}",
            help=msg,
        )

    parser.add_argument("--prefix", "--pf", type=str, default=None, help="file prefix")

    parser.add_argument("--output", "--o", type=str, default=None, help="Quiet mode")

    # Column selection
    parser.add_argument(
        "--columns",
        type=str,
        default=None,
        nargs="*",
        help="Specify columns to read (e.g., '1,2,5' or 'colA,colB' or '1 2' or 'Time Lx,Lz' etc.)",
    )
    parser.add_argument(
        "--select",
        "--sel",
        type=str,
        nargs="+",
        action="extend",
        help="Range selection for a column",
    )
    parser.add_argument(
        "--skim",
        "--skip",
        type=int,
        default=None,
        help="Skimming data interval",
    )
    parser.add_argument(
        "--plottype",
        "--pt",
        type=str,
        default=None,
        choices=["plotly", "pl", "matplotlib", "pyplot", "pp"],
        help="Range selection for a column",
    )

    return parser


def combine_inputs(args):
    """
    Combine all input_* attributes from argparse Namespace into args.input.

    For each input_* attribute:
    - If it has a suffix (e.g., input_txt, input_dat), use that as the type
    - If it's just input_ with no suffix, extract the file extension
    - Creates a list of dicts with 'filename' and 'type' keys

    Args:
        args: argparse Namespace object
    """
    combined = []

    # Iterate through all attributes in the namespace
    for attr_name in vars(args):
        if not attr_name.startswith("input"):
            continue

        values = getattr(args, attr_name)

        # Skip if None or empty
        if not values:
            continue

        # Ensure values is a list
        if not isinstance(values, list):
            values = [values]

        # Determine the type
        if attr_name == "input_":
            # No suffix, extract extension from filename
            for filename in values:
                ext = Path(filename).suffix.lstrip(".")
                combined.append(
                    {
                        "filename": filename,
                        "type": ext,
                    }
                )
        else:
            # Has suffix (e.g., input_txt, input_dat, input_csv)
            file_type = attr_name.split("_", 1)[1]  # Get 'txt', 'dat', 'csv', etc.
            for filename in values:
                combined.append(
                    {
                        "filename": filename,
                        "type": file_type,
                    }
                )

    # Assign to args.input
    args.input = combined
    return combined
