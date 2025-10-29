#!/usr/bin/env python3
"""Execute shell commands with optional progress tracking and logging."""

import os
import re
import sys
import logging
import subprocess
from typing import List, Union, Optional


def execute_commands(
    commands: Union[str, List[str]],
    dry_run: bool = False,
    logger: Optional[Union[logging.Logger, str]] = None,
) -> List[Union[subprocess.CompletedProcess, str]]:
    """
    Execute a list of shell commands and return the results.

    Parameters
    ----------
    commands : str or list of str or list of list of str
        Shell command(s) to execute. Can be:
        - A single command string
        - A list of command strings
        - A list of lists of command strings (nested structure)
    dry_run : bool, optional
        If True, log commands without executing them. Default is False.
    logger : logging.Logger or str, optional
        Logger instance or logging level string ('DEBUG', 'INFO', etc.).
        If None, defaults to INFO level.

    Returns
    -------
    list
        List of subprocess.CompletedProcess objects for each command.
        Returns list of "Dry run" strings if dry_run is True.

    Raises
    ------
    ValueError
        If commands is not a string or list of strings.
    SystemExit
        If any command fails (non-zero return code).
    """
    # Setup logging
    if isinstance(logger, str):
        level = getattr(logging, logger.upper(), logging.INFO)
    elif logger is None:
        level = logging.INFO
    else:
        level = logger.level if hasattr(logger, "level") else logging.INFO

    logging.basicConfig(level=level, format="%(levelname)s - %(message)s")
    log = logger if isinstance(logger, logging.Logger) else logging.getLogger(__name__)

    # Validate and normalize commands input
    if isinstance(commands, str):
        commands = [commands]
    elif not isinstance(commands, list):
        raise ValueError("Commands must be a string or list of strings")

    # Handle dry run
    if dry_run:
        try:
            for command in commands:
                log.info(f"  [DRY RUN] {re.sub(r'\s+', ' ', command)}")
        except TypeError:
            for command_list in commands:
                for command in command_list:
                    log.info(f"  [DRY RUN] {re.sub(r'\s+', ' ', command)}")

        return ["Dry run"] * len(commands)

    # Execute commands
    results = []
    with _get_progress_bar(len(commands)) as bar:
        try:
            # Try to execute as flat list of commands
            for command in commands:
                result = _execute_single_command(command, log)
                results.append(result)
                bar()
        except TypeError:
            # Handle nested list structure (list of lists)
            results = []
            for command_list in commands:
                for command in command_list:
                    result = _execute_single_command(command, log)
                    results.append(result)
                bar()

    return results


def _execute_single_command(
    command: str, logger: logging.Logger
) -> Union[subprocess.CompletedProcess, str]:
    """
    Execute a single shell command.

    Parameters
    ----------
    command : str
        Shell command to execute.
    logger : logging.Logger
        Logger instance for output.

    Returns
    -------
    subprocess.CompletedProcess or str
        Result of command execution, or empty string for cd commands.
    """
    logger.info(f"  Running: {re.sub(r'\s+', ' ', command)}")

    # Handle cd command specially
    if command.startswith("cd "):
        target_dir = command[3:].strip()
        os.chdir(target_dir)
        logger.debug(f"  Changed directory to: {os.getcwd()}")
        logger.debug("-" * 80)
        return ""

    # Execute regular command
    result = subprocess.run(
        command, shell=True, text=True, capture_output=True, check=False
    )

    if result.returncode != 0:
        logger.error("  Command failed")
        for line in result.stderr.strip().split("\n"):
            if line:  # Skip empty lines
                logger.error(f"  {line}")
        sys.exit(1)
    else:
        for line in result.stdout.strip().split("\n"):
            if line:  # Skip empty lines
                logger.debug(f"  {line}")

    logger.debug("-" * 80)
    return result


def _get_progress_bar(total: Optional[int] = None):
    """
    Get a progress bar context manager.

    Uses alive_progress if available, otherwise returns a dummy implementation.

    Parameters
    ----------
    total : int, optional
        Total number of items for progress tracking.

    Returns
    -------
    Context manager
        Progress bar context manager (or dummy if alive_progress not available).
    """
    try:
        from alive_progress import alive_bar, config_handler

        config_handler.set_global(
            force_tty=True, receipt=True, monitor=True, bar="smooth"
        )
        return alive_bar(total)
    except ImportError:
        return _DummyProgressBar(total)


class _DummyProgressBar:
    """Fallback progress bar when alive_progress is not available."""

    def __init__(self, total: Optional[int] = None):
        self.total = total

    def __enter__(self):
        return self.__call__

    def __exit__(self, *args):
        pass

    def __call__(self, *args, **kwargs):
        pass
