#!/usr/bin/env python

"""
Quick plotting utilities for scientific data visualization.

This module provides utilities for reading various file formats and
creating publication-quality plots with matplotlib.
"""

import sys
import pandas as pd
import logging

from .logger import setup_logger
from .basic_parser import basic_parser_2, combine_inputs
from .read_files import read_multi_files_to_df
from .plotting_wrapper import plot
from .fitting import fitting_wrapper


def main():
    """Main function for command-line quick plotting."""
    parser = basic_parser_2(description="Quick plotting utility")

    parser.add_argument("--xlabel", "--xl", type=str, default=None, help="x-axis label")
    parser.add_argument("--ylabel", "--yl", type=str, default=None, help="y-axis label")

    parser.add_argument(
        "--logx", "--lx", action="store_true", help="log scale for x-axis"
    )
    parser.add_argument(
        "--logy", "--ly", action="store_true", help="log scale for x-axis"
    )

    parser.add_argument("--legend", nargs="*", help="Legend labels for the data series")

    parser.add_argument("--fit", type=str, default=None, help="Add a fit to the data")

    # Parse arguments
    args, unknown = parser.parse_known_args()

    # Set up logging
    logger = setup_logger(args=args)
    logger.info("Quick plotting program")

    if unknown:
        args = parse_unknown_as_gnuplot(unknown, args)

    combine_inputs(args)
    if len(args.input) == 0:
        logger.error("No input files were specified")
        sys.exit(1)
    logger.debug(args.input)

    # Process columns argument, allow for gnuplot-style colon separator
    if args.columns is None:
        args.columns = ["1", "2"]
    else:
        if len(args.columns) == 1 and ":" in args.columns[0]:
            args.columns = args.columns[0].split(":")
    logger.debug(f"... Selected columns: {args.columns}")

    # Read data
    columns_selected, df_selected = read_multi_files_to_df(
        args.input,
        mask=args.select,
        columns=args.columns,
        prefix=args.prefix,
        skim=args.skim,
    )

    # Separete x and y data for processing
    x_data = df_selected[0]
    y_data = df_selected[1]

    # Defile axis labels
    # xlabel = "Index"
    # ylabel = "Values"
    if len(columns_selected) == 1:
        ylabel = columns_selected[0][0] if len(columns_selected[0]) == 1 else "Y values"

    elif len(columns_selected) == 2:
        xlabel = columns_selected[0][0] if len(columns_selected[0]) == 1 else "X values"
        ylabel = columns_selected[1][0] if len(columns_selected[1]) == 1 else "Y values"
    else:
        logger.error(
            f"Quick plot - don't know what to do with a {len(columns_selected)} dimensional data set"
        )
    xlabel = args.xlabel if args.xlabel is not None else xlabel
    ylabel = args.ylabel if args.ylabel is not None else ylabel
    logger.debug(f"... x-axis label: {xlabel}")
    logger.debug(f"... y-axis label: {ylabel}")

    # Fit data and add fitted curves
    if args.fit is not None:
        res = fitting_wrapper(x_data, y_data, fit_type=args.fit)

        fitted_x_data = []
        fitted_y_data = []

        for i in range(len(y_data.columns)):
            # Create unique column names for fitted x data
            fitted_x_col_name = (
                f"fitted_range_{i}" if len(y_data.columns) > 1 else "fitted_range"
            )

            string = res[i]["evaluated_expression"].replace("+ -", "-")
            string = string.replace("*x", "x")
            fitted_x_data.append(
                pd.DataFrame(
                    res[i]["fitted_range"],
                    columns=[fitted_x_col_name],
                )
            )

            fitted_y_data.append(
                pd.DataFrame(
                    res[i]["fitted_values"],
                    columns=[string],
                )
            )
            columns_selected[0].append(fitted_x_col_name + "_fit")
            columns_selected[1].append(string)

        # Concatenate all fitted data at once
        x_data = pd.concat([x_data] + fitted_x_data, axis=1)
        y_data = pd.concat([y_data] + fitted_y_data, axis=1)

    # Add legend if multiple data series
    legend = args.legend
    if legend is None:
        if len(columns_selected[1]) > 1:
            legend = [f"{c}" for c in columns_selected[1]]
        elif len(args.input) > 1:
            legend = [f"{c['filename']}" for c in args.input]
    logger.debug(f"Legend: {legend}")

    # Set plot title
    if len(args.input) > 1:
        ptitle = "Plot of multiple files"
        if args.prefix:
            ptitle += f" in {args.prefix}"
    else:
        prefix = args.prefix + "/" if args.prefix else ""
        fname = f"{prefix}{args.input[0]['filename']}".replace("//", "/")
        if len(fname) > 50:
            fname = "...." + fname[-47:]
        ptitle = f"Plot of {fname}"

    # Select plot type and file output for matplotlib
    fout = args.output
    if args.plottype is None:
        ptype = "plotly"
    elif args.plottype in ["pl", "plotly"]:
        ptype = "plotly"
    elif args.plottype in ["pp", "pyplot", "matplotlib"]:
        ptype = "matplotlib"
    logger.debug(f"... Plot type {ptype}")
    if fout:
        logger.debug(f"... Output filename {fout}")

    # Create the plot
    plot(
        x_data,
        y_data,
        xlabel=xlabel,
        ylabel=ylabel,
        title=ptitle,
        legend=legend,
        ptype=ptype,
        fout=fout,
        grid=True,
        logx=args.logx,
        logy=args.logy,
    )


def parse_unknown_as_gnuplot(unknown, args):
    logger = logging.getLogger("mylogger")
    logger.warning("Parsing unknow commands as gnuplot. Mixed input may not work")
    logger.debug(f"Unknown string: {' '.join(unknown)}")

    args.input_ = [unknown[0]]
    logger.debug(f"Input file assumed to be: {unknown[0]}")
    if "u" in unknown:
        args.columns = [unknown[unknown.index("u") + 1]]
        logger.debug(f"Columns to plot assumed to be: {args.columns[0]}")

    return args


if __name__ == "__main__":
    main()
