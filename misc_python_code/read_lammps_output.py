#!/usr/bin/env python3
"""
LAMMPS output file parser.

Reads LAMMPS output files and extracts data blocks into pandas DataFrames.
"""
import io
import sys
import logging
from pathlib import Path
from typing import List, Optional

import pandas as pd


def read_lammps_output_to_df(
    file_path: str,
    logger: Optional[logging.Logger] = None,
    header: str = "Step",
    footer: str = "Loop",
) -> List[pd.DataFrame]:
    """
    Read LAMMPS output file and extract data blocks into DataFrames.
    
    Args:
        file_path: Path to the LAMMPS output file
        logger: Logger instance (creates one if None)
        header: String that marks the start of a data block
        footer: String that marks the end of a data block
    
    Returns:
        List of DataFrames containing the extracted data blocks
    
    Raises:
        FileNotFoundError: If the input file doesn't exist
        ValueError: If no valid data blocks are found
        PermissionError: If the file cannot be read due to permissions
    """
    # Set up logger if not provided
    if logger is None:
        logging.basicConfig(
            level=logging.INFO,
            format='%(levelname)s - %(message)s'
        )
        logger = logging.getLogger(__name__)
    
    # Validate file path
    file_path = Path(file_path)
    if not file_path.exists():
        error_msg = f"File not found: {file_path}"
        logger.error(error_msg)
        raise FileNotFoundError(error_msg)
    
    if not file_path.is_file():
        error_msg = f"Path is not a file: {file_path}"
        logger.error(error_msg)
        raise ValueError(error_msg)
    
    logger.info(f"Reading LAMMPS output file: {file_path}")
    
    # Read file with error handling
    try:
        with open(file_path, "r", encoding='utf-8') as file:
            lines = file.readlines()
    except PermissionError as e:
        error_msg = f"Permission denied reading file: {file_path}"
        logger.error(error_msg)
        raise PermissionError(error_msg) from e
    except UnicodeDecodeError:
        logger.warning("UTF-8 decoding failed, trying with latin-1 encoding")
        try:
            with open(file_path, "r", encoding='latin-1') as file:
                lines = file.readlines()
        except Exception as e:
            error_msg = f"Failed to read file with multiple encodings: {e}"
            logger.error(error_msg)
            raise ValueError(error_msg) from e
    except Exception as e:
        error_msg = f"Unexpected error reading file: {e}"
        logger.error(error_msg)
        raise IOError(error_msg) from e
    
    logger.debug(f"Read {len(lines)} lines from file")
    
    if not lines:
        error_msg = "File is empty"
        logger.error(error_msg)
        raise ValueError(error_msg)
    
    # Extract data blocks
    data_blocks = []
    i = 0
    block_count = 0
    
    while i < len(lines):
        line_stripped = lines[i].strip()
        
        if line_stripped.startswith(header):
            block_count += 1
            header_index = i
            logger.debug(f"Found '{header}' header at line {i} (block {block_count})")
            
            # Find the end of this block
            footer_index = None
            for j in range(header_index + 1, len(lines)):
                line_j_stripped = lines[j].strip()
                if line_j_stripped.startswith(footer) or line_j_stripped.startswith("WARNING"):
                    footer_index = j
                    logger.debug(f"Block {block_count} ends at line {footer_index}")
                    break
            
            # Extract data lines for this block
            if footer_index is None:
                data_lines = lines[header_index:]
                i = len(lines)
                logger.debug(f"Block {block_count} extends to end of file")
            else:
                data_lines = lines[header_index:footer_index]
                i = footer_index + 1
            
            # Parse block into DataFrame
            data_text = "".join(data_lines)
            try:
                df = pd.read_csv(io.StringIO(data_text), sep=r"\s+")
                
                if df.empty:
                    logger.warning(f"Block {block_count} is empty, skipping")
                else:
                    data_blocks.append(df)
                    logger.debug(
                        f"Successfully parsed block {block_count} with shape {df.shape}"
                    )
            except pd.errors.EmptyDataError:
                logger.warning(f"Block {block_count} contains no data, skipping")
            except pd.errors.ParserError as e:
                logger.warning(f"Failed to parse data block {block_count}: {e}")
            except Exception as e:
                logger.warning(
                    f"Unexpected error parsing block {block_count}: {type(e).__name__}: {e}"
                )
        else:
            i += 1
    
    # Validate results
    if not data_blocks:
        error_msg = f"Could not find any valid data blocks starting with '{header}'"
        logger.error(error_msg)
        raise ValueError(error_msg)
    
    logger.info(f"Successfully extracted {len(data_blocks)} data block(s)")
    return data_blocks


def print_dataframe_info(df: pd.DataFrame, block_num: int) -> None:
    """
    Print information about a DataFrame's structure.
    
    Args:
        df: DataFrame to inspect
        block_num: Block number for labeling
    """
    print(f"\n{'=' * 70}")
    print(f"Block {block_num} Structure")
    print(f"{'=' * 70}")
    print(f"Shape: {df.shape[0]} rows × {df.shape[1]} columns")
    print(f"\nColumns: {', '.join(df.columns.tolist())}")
    print(f"\nData types:")
    print(df.dtypes)
    print(f"\nFirst few rows:")
    print(df.head())
    print(f"\nBasic statistics:")
    print(df.describe())


def main():
    """Main entry point for the script."""
    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    logger = logging.getLogger(__name__)
    
    # Check command line arguments
    if len(sys.argv) < 2:
        print("Usage: python script.py <lammps_output_file>", file=sys.stderr)
        print("\nExample: python script.py output.lammps", file=sys.stderr)
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    try:
        # Read LAMMPS output file
        data_blocks = read_lammps_output_to_df(input_file, logger=logger)
        
        # Display structure of each data block
        print(f"\nFound {len(data_blocks)} data block(s) in file: {input_file}")
        
        for i, df in enumerate(data_blocks, start=1):
            print_dataframe_info(df, i)
        
        # Summary
        print(f"\n{'=' * 70}")
        print("Summary")
        print(f"{'=' * 70}")
        print(f"Total blocks: {len(data_blocks)}")
        print(f"Total rows: {sum(df.shape[0] for df in data_blocks)}")
        
    except FileNotFoundError as e:
        logger.error(f"File error: {e}")
        sys.exit(1)
    except ValueError as e:
        logger.error(f"Value error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {type(e).__name__}: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
