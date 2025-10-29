#!/usr/bin/env python

from ..logger import setup_logger, formatting
from .commad_line_parser import setup_argument_parser
from .analyse_orca import process_orca_file


def main():
    args = setup_argument_parser()

    logger = setup_logger(args=args)

    logger.debug("-" * 50)
    for k, v in args.__dict__.items():
        logger.debug(f"{k}: {v}")
    logger.info(formatting().dashes)

    process_orca_file(args.__dict__)


if __name__ == "__main__":
    main()
