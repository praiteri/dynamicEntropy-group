"""
Quick plotting utilities for scientific data visualization.

This package provides utilities for reading various file formats and
creating publication-quality plots with matplotlib and plotly.
"""

from .logger import *
from .basic_parser import *
from .read_files import *

from .plotting_wrapper import *
from .plotting_utilities import *
from .plotting_with_matplotlib import *
from .plotting_with_plotly import *

from .plot import *

from .fitting import *

__version__ = "0.1.0"
__author__ = "Paolo Raiteri"
__email__ = "p.raiteri@curtin.edu.au"

__all__ = [
    "setup_logger",
    "test_logger",
    "read_file_to_df",
    "plot",
]
