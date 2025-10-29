#!/usr/bin/env python

"""
Quick plotting utilities for scientific data visualization.

This module provides utilities for reading various file formats and
creating publication-quality plots with matplotlib.
"""

import sys

import numpy as np
import pandas as pd

import logging
from .logger import setup_logger
from .basic_parser import basic_parser, basic_parser_2, combine_inputs
from .read_files import read_file_to_df, read_multi_files_to_df
from .select_columns import select_columns
from .plotting_wrapper import plot


def normalized_gaussian(x, mu, sigma):
    """
    Calculates the normalized Gaussian (Normal Distribution PDF).

    Args:
        x (numpy.ndarray or float): The input value(s) at which to evaluate the PDF.
        mu (float): The mean of the Gaussian distribution.
        sigma (float): The standard deviation of the Gaussian distribution.

    Returns:
        numpy.ndarray or float: The probability density at the given x value(s).
    """
    coefficient = 1 / (sigma * np.sqrt(2 * np.pi))
    exponent = -((x - mu) ** 2) / (2 * sigma**2)
    return coefficient * np.exp(exponent)


def compute_histogram(data, bins="auto", range=None, weights=None, density=False):
    """
    Compute histogram of a dataset using NumPy.

    Parameters:
    -----------
    data : array-like
        Input data to compute histogram
    bins : int or str, optional (default='auto')
        Number of bins or method to compute optimal bins
        Can be an integer or one of ['auto', 'fd', 'scott', 'rice', 'sturges', 'sqrt']
    range : tuple, optional (default=None)
        The lower and upper range of the bins
    weights : array-like, optional (default=None)
        Weights for each data point
    density : bool, optional (default=False)
        If True, return probability density function

    Returns:
    --------
    counts : array
        The counts or density for each bin
    bin_edges : array
        The edges of each bin
    bin_centers : array
        The centers of each bin
    """
    # Compute histogram
    counts, bin_edges = np.histogram(data, bins=bins, range=range, weights=weights, density=density)

    # Calculate bin centers for easier plotting
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    return counts, bin_edges, bin_centers


def main():
    """Main function for command-line quick plotting."""
    parser = basic_parser_2(description="Quick plotting utility")

    # Parse arguments
    args, unknown = parser.parse_known_args()
    combine_inputs(args)

    # Set up logging
    logger = setup_logger(args=args)
    logger.info("Quick histogram program")

    # Select plot type and file output for matplotlib
    fout = None
    if hasattr(args, "output"):
        fout = args.output if args.output else None

    ptype = "plotly"
    if hasattr(args, "plottype"):
        if args.plottype in ["pp", "pyplot", "matplotlib"]:
            ptype = "matplotlib"

    # # Read data
    # data = read_file_to_df(args.input, mask=args.select)

    # # Parse column indices or names
    # _, df_selected = select_columns(args, data)
    columns_selected, df_selected = read_multi_files_to_df(
        args.input, mask=args.select, columns=args.columns, prefix=args.prefix
    )

    if len(df_selected) > 1:
        df = pd.concat(df_selected, axis=1)
    else:
        df = df_selected[0]

    xlabel = "Values"
    ylabel = "Probability"
    x_data = []
    y_data = []
    customstyles = []
    legend = []
    for col_name, col_data in df.items():
        histo, bin_edges, bin_centres = compute_histogram(col_data)
        x_data.append(bin_centres)
        y_data.append(histo)
        legend.append(col_name)
        customstyles.append("b")

    plot(
        x_data,
        y_data,
        xlabel=xlabel,
        ylabel=ylabel,
        title=f"Histogram",
        legend=legend,
        ptype=ptype,
        fout=fout,
        grid=True,
        logx=False,
        logy=False,
        customstyles=customstyles,
    )


if __name__ == "__main__":
    main()
