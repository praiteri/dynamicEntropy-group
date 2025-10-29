#!/usr/bin/env python

"""
Quick plotting utilities for scientific data visualization.

This module provides utilities for reading various file formats and
creating publication-quality plots with matplotlib.
"""
import os
import sys
import argparse
import logging

from .logger import setup_logger

from .basic_stats import main as stats
from .plot import main as plot
from .histogram import main as histogram
from .kb_integral import main as kb
from .association_constant import main as pkdiss


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--wrap",
        "--w",
        type=str,
        default=None,
        nargs="*",
        help="Wrapper to data processing programs",
    )

    args, unknown = parser.parse_known_args()

    if args.wrap is None:
        print("Nothing to do")
        sys.exit(0)

    cmds = [item for x in args.wrap for item in x.split(",")]

    for c in cmds:
        print(c)
        if c.lower() in "statistics":
            stats()
        elif c.lower() in "plotting":
            plot()
        elif c.lower() in "histogram":
            histogram()
        elif c.lower() in "kb_integral":
            kb()
        elif c.lower() in "pkdiss":
            pkdiss()
        else:
            print(f"Uknown program: {c}")
            sys.exit(1)

    # logger = setup_logger()
    # logger.info("Data processing wrapper")

    # if unknown[0] in "histogram":
    #     histogram()


if __name__ == "__main__":
    main()
