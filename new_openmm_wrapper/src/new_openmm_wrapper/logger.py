### Copyright (C) 2026  Paolo Raiteri
### This program is free software: you can redistribute it and/or modify
### it under the terms of the GNU General Public License as published by
### the Free Software Foundation, either version 3 of the License, or
### (at your option) any later version.

import os
import sys
import logging
import colorama
from colorama import Fore, Style

# Initialize colorama to work on Windows too
colorama.init()


class ColoredFormatter(logging.Formatter):
    """Custom formatter that adds colors to log messages based on level with fixed width"""

    COLORS = {
        "DEBUG": Fore.CYAN + Style.BRIGHT,
        "INFO": Fore.WHITE,
        "WARNING": Fore.YELLOW + Style.BRIGHT,
        "ERROR": Fore.RED + Style.BRIGHT,
        "CRITICAL": Fore.MAGENTA + Style.BRIGHT,
    }
    LEVEL_NAME_WIDTH = 8

    def format(self, record):
        # Save original levelname and message
        orig_levelname = record.levelname
        orig_msg = record.msg

        # Determine the color for the level
        color_code = self.COLORS.get(orig_levelname, "")

        # Calculate padding for levelname
        padding = " " * (self.LEVEL_NAME_WIDTH - len(orig_levelname))

        # Apply color to levelname with padding
        record.levelname = f"{color_code}{orig_levelname}{padding}{Style.RESET_ALL}"

        # Apply color to the message as well
        record.msg = f"{color_code}{orig_msg}{Style.RESET_ALL}"

        # Format the message
        result = super().format(record)

        # Restore original values
        record.levelname = orig_levelname
        record.msg = orig_msg
        return result


def create_logger(level="INFO", filename=None):
    """Set up and return a colored logger instance"""

    # Create formatter
    formatter = ColoredFormatter(fmt="%(levelname)8s - %(message)s")

    # Create logger
    logger = logging.getLogger("dynamicEntropy")

    # Set logger level
    try:
        logger.setLevel(getattr(logging, level.upper()))
    except AttributeError:
        logger.setLevel(logging.INFO)  # Default to INFO if invalid level

    # Create console handler
    if filename is None:
        handler = logging.StreamHandler(sys.stdout)

        # Add formatter to handler
        handler.setFormatter(formatter)

    else:
        if os.path.exists(filename):
            os.remove(filename)
        handler = logging.FileHandler(filename)

    # Set handler level
    handler.setLevel(level)

    # Add handler to logger
    logger.addHandler(handler)

    return logger
