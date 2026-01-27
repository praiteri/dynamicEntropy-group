### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import logging

indentation = 2


def pretty_log(
    data=None,
    logger=None,
    indent=0,
    title=None,
    marker=":",
    align_width=None,
    sep=False,
    ncols=6,
):
    if logger == "print":
        logger = print
    else:
        if logger is None:
            logger = "INFO"

        logger_obj = logging.getLogger("dynamicEntropy")
        if logger.upper() == "DEBUG":
            logger = logger_obj.debug
        elif logger.upper() == "INFO":
            logger = logger_obj.info
        elif logger.upper() == "WARNING":
            logger = logger_obj.warning
        elif logger.upper() == "ERROR":
            logger = logger_obj.error
        elif logger.upper() == "CRITICAL":
            logger = logger_obj.critical
        else:
            logger = logger_obj.info

    if sep:
        logger("=" * 50)

    if title is not None:
        indent_str = " " * (indent * indentation)
        logger(f"{indent_str}{title}")

    if data is None:
        return

    elif isinstance(data, dict):
        if title is not None:
            indent += 1
        _pretty_log_dict(data, logger, indent, marker, align_width, ncols)

    elif isinstance(data, list):
        if title is not None:
            indent += 1
        _pretty_log_list(data, logger, indent, ncols)

    else:
        indent_str = " " * (indent * indentation)
        if title is None:
            logger(f"{indent_str}{data}")
        else:
            indent_str = " " * (indent + 2)
            logger(f"{indent_str}{title:<}{marker} {data}")


def _pretty_log_list(data, logger, indent, ncols=6):
    indent_str = " " * (indent * indentation)
    n = len(data)
    for i in range(0, n, ncols):
        string = " ".join(f"{x:8}" for x in data[i : i + ncols])
        logger(f"{indent_str}{string}")


def _pretty_log_dict(
    data, logger=print, indent=0, marker=":", align_width=None, ncols=6
):
    """
    Format a dictionary as a pretty-printed string with indentation.

    Each nested level is indented by 2 spaces.

    Args:
        data: Dictionary or value to format
        logger: Function to call for output (default: print)
        indent (int): Current indentation level (number of spaces)
    """
    indent_str = " " * (indent * indentation)

    if align_width is None:
        align_col = 48 - len(indent_str)
    else:
        align_col = max(0, align_width - len(indent_str))

    for key, value in data.items():
        if isinstance(value, dict):
            # Nested dictionary: write key and recurse
            logger(f"{indent_str}{key:<}{marker}")
            _pretty_log_dict(
                value, logger, indent + 1, ncols=ncols, align_width=align_width
            )
        elif isinstance(value, (str, int, float)):
            # String: write as-is
            logger(f"{indent_str}{key:<{align_col}} {marker} {value}")
        elif isinstance(value, bool):
            # Boolean: write as lowercase true/false (YAML style)
            logger(f"{indent_str}{key:<{align_col}} {marker} {str(value).lower()}")
        elif value is None:
            # None: write as null
            logger(f"{indent_str}{key:<{align_col}} {marker} null")
        elif isinstance(value, list):
            # List: write inline
            logger(f"{indent_str}{key}{marker}")
            _pretty_log_list(value, logger, indent + 1, ncols=ncols)
        else:
            # Numbers and other types: write directly
            logger(f"{indent_str}{key:<{align_col}} {marker} {value}")

    return


def test_pretty_log():
    test_dict = {
        "name": "water box",
        "num_atoms": 512,
        "temperature": 300.0,
    }
    test_list = [1, 2, 3, 4, 5]

    # Separator
    pretty_log(logger="debug")
    pretty_log()
    pretty_log(logger="info")
    pretty_log(logger="warning")
    pretty_log(logger="error")
    pretty_log(logger="critical")
    print("***")
    # Text only
    pretty_log(None, title="charge", logger="INFO", marker="->")
    print("***")
    # variable and value
    pretty_log(12, title="charge", logger="INFO")
    pretty_log(12, title="charge", logger="INFO", marker=":")
    pretty_log(12, title="charge", logger="INFO", marker="->")
    print("***")
    # variable and value as dict
    pretty_log({"charge": 12}, logger="INFO", marker=":")
    print("***")
    # Title and dict
    pretty_log(test_dict, logger="INFO", title="test dictionary:")
    pretty_log(
        test_dict,
        logger="DEBUG",
    )
    print("***")
    pretty_log(test_list * 5, title="test list:")

    quit()
