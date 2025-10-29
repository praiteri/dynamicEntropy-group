#!/usr/bin/env python

"""
Quick plotting utilities for scientific data visualization.

This module provides utilities for reading various file formats and
creating publication-quality plots with matplotlib.
"""

import numpy as np
import pandas as pd
import scipy.stats as stats
import logging

from .logger import setup_logger
from .basic_parser import basic_parser
from .read_files import read_file_to_df
from .select_columns import select_columns
from .format_output_dict import format_output_dict_auto


def basic_stats(df, **kwargs):
    logger = logging.getLogger("mylogger")

    logger.info("Computing basic statistics ...")
    for c in df.columns:
        stats = compute_statistics(np.array(df[c]))

        if "minimal" in kwargs:
            if type(c) is int:
                msg = str(c)
            else:
                msg = c
            avg = stats["Average"]
            ci = stats["95% Confidence Interval"]
            kr = stats["Kurtosis"]
            logger.result(f" {msg:>30s} : {avg:>13.5g} +/- {ci:<13.5g} | {kr:>6.3f}")
        else:
            msg = "Column name"
            logger.result(f"{msg:>25s} : {c}")
            output = format_output_dict_auto(stats)
            for l in output:
                logger.result(l)
            # for k, v in stats.items():
            #     logger.result(f"{k:>25s} : {v}")
            if len(df.columns) > 1:
                logger.result("-" * 50)


def compute_statistics(data):
    # Critical value from the t-distribution
    confidence_level = 0.95
    t_critical = stats.t.ppf((1 + confidence_level) / 2, (len(data) - 1))

    def drift(data):
        xx = np.array([i for i in range(len(data))])
        linear_fit = np.polyfit(xx, data, deg=1)
        return linear_fit[0]

    statistics = {
        "Number of values": len(data),
        "Average": np.average(data, axis=0),
        "Standard Deviation": np.std(data, axis=0),
        "95% Confidence Interval": t_critical * np.std(data, axis=0) / np.sqrt(len(data)),
        "Minimum value": np.min(data, axis=0),
        "Maximum value": np.max(data, axis=0),
        "Skewness": stats.skew(data),
        "Kurtosis": stats.kurtosis(data),
        "Drift (last-first)": (data[-1] - data[0]),
        "Drift (linear fit)": drift(data) * len(data),
    }

    return statistics


def main():
    """Main function for command-line quick statistics."""
    parser = basic_parser(description="Quick statistics utility")

    # Parse arguments
    args, unknown = parser.parse_known_args()

    # Set up logging
    logger = setup_logger(args=args)
    logger.info("Computing basic statistics")

    # Read data
    data = read_file_to_df(args.input, mask=args.select)

    # Parse column indices or names
    _, df_selected = select_columns(args.columns, data)

    if len(df_selected) > 1:
        df = pd.concat(df_selected, axis=1)
    else:
        df = df_selected[0]

    basic_stats(df)


if __name__ == "__main__":
    main()
