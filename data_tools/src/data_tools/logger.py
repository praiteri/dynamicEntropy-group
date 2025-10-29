import os
import sys
import logging
import colorama
from colorama import Fore, Style

# ============================================================================
# Logger Setup
# ============================================================================


class formatting:
    def __init__(self):
        self.ldash = 120
        self.ndash = 80
        self.sdash = 40
        self.dashes = "-" * self.ndash
        self.dashes_short = "-" * self.sdash
        self.dashes_long = "-" * self.ldash
        self.fmt = "{:>40} :: {:<36}"


class ColoredFormatter(logging.Formatter):
    """Custom formatter that adds colors to log messages based on level with fixed width"""

    COLORS = {
        "DEBUG": Style.BRIGHT + Fore.CYAN,
        "VERBOSE": Style.BRIGHT + Fore.BLUE,
        "INFO": Style.BRIGHT + Fore.WHITE,
        "RESULT": Style.BRIGHT + Fore.GREEN,
        "WARNING": Style.BRIGHT + Fore.YELLOW,
        "ERROR": Style.BRIGHT + Fore.RED,
        "CRITICAL": Style.BRIGHT + Fore.MAGENTA,
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


def setup_logger(logger_name="mylogger", level="INFO", filename=None, args=None):
    """Set up and return a colored logger instance"""

    # Avoid creating a logger if it aready exists
    if logger_name in logging.Logger.manager.loggerDict:
        return logging.getLogger(logger_name)

    # Define custom VERBOSE level (between DEBUG-10 and INFO-20)
    VERBOSE = 15
    RESULT = 25

    class ColoredLogger(logging.Logger):
        """Custom logger class with verbose method"""

        def verbose(self, msg, *args, **kwargs):
            """Log at custom VERBOSE level"""
            # if self.isEnabledFor(VERBOSE):
            self._log(VERBOSE, msg, args, **kwargs)

        def result(self, msg, *args, **kwargs):
            """Log at custom RESULT level"""
            # if self.isEnabledFor(RESULT):
            self._log(RESULT, msg, args, **kwargs)

    # Initialize colorama to work on Windows too
    colorama.init()

    # Create formatter
    formatter = ColoredFormatter(
        # fmt='%(asctime)s - %(levelname)s - %(message)s',
        # datefmt='%Y-%m-%d %H:%M:%S'
        fmt="%(levelname)8s - %(message)s"
    )

    # Register custom levels
    logging.addLevelName(VERBOSE, "VERBOSE")
    logging.addLevelName(RESULT, "RESULT")
    logging.VERBOSE = VERBOSE
    logging.RESULT = RESULT

    # Register our custom logger class
    logging.setLoggerClass(ColoredLogger)

    # Create logger
    logger = logging.getLogger(logger_name)

    # Use argparse Namespace for settings
    if args is not None:
        # Process log levels (first valid option wins)
        if hasattr(args, "debug") and args.debug:
            level = "DEBUG"
        elif hasattr(args, "verbose") and args.verbose:
            level = "VERBOSE"
        elif hasattr(args, "result") and args.quiet:
            level = "RESULT"
        elif hasattr(args, "quiet") and args.quiet:
            level = "ERROR"

        # Process log filename
        if (log_path := getattr(args, "logfile", None)) is not None:
            filename = log_path

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

    # To avoid with libraries that configure logging without checking if it's already set up.
    logger.propagate = False

    # Set handler level
    handler.setLevel(level)

    # Add handler to logger
    logger.addHandler(handler)

    # test_logger()

    return logger


def test_logger():
    logger = logging.getLogger("mylogger")
    logger.debug("This is a debug message")
    logger.verbose("This is an verbose message")
    logger.info("This is an info message")
    logger.result("This is an result message")
    logger.warning("This is a warning message")
    logger.error("This is an error message")
    logger.critical("This is a critical message")
