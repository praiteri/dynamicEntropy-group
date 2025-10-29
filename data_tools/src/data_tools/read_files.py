import os
import io
import sys

from typing import Optional, Union, List, Dict, Tuple

import tempfile
import pandas as pd
import subprocess as sp

import logging
from .logger import formatting
from .mask_df_by_range import mask_df_by_range
from .select_columns import select_columns


# ============================================================================
def is_numeric(fields):
    for x in fields:
        x = x.lower().strip()
        # Handle inf and nan explicitly
        if x not in ("inf", "-inf", "nan", "+inf", "infinity", "-infinity"):
            try:
                float(x)
            except ValueError:
                return False
        return True


def skim_dataframe(df, n):
    """
    Return every nth row from a DataFrame.

    Parameters:
    -----------
    df : pandas.DataFrame
        The input DataFrame to skim
    n : int
        The step size - returns every nth row

    Returns:
    --------
    pandas.DataFrame
        A new DataFrame containing every nth row

    Examples:
    ---------
    >>> import pandas as pd
    >>> df = pd.DataFrame({'A': range(10), 'B': range(10, 20)})
    >>> skim_dataframe(df, 3)  # Returns rows 0, 3, 6, 9
    """
    return df.iloc[::n]


# ============================================================================
# File Reading Utilities
# ============================================================================


def read_multi_files_to_df(
    list_of_files: Union[List[Union[str, Dict[str, str]]], Union[str, Dict[str, str]]],
    mask: Optional[Dict[str, Tuple[float, float]]] = None,
    columns: Optional[List[Union[int, str]]] = None,
    skim: Optional[int] = None,
    prefix: Optional[str] = "",
) -> Tuple[Optional[List[Union[int, str]]], List[pd.DataFrame]]:
    """Read one or more files into a pandas DataFrame.

    Parameters
    ----------
    list_of_files : str, dict, or list of str/dict
        Path(s) to the file(s), or dict(s) with 'filename' and 'type' keys
    mask : dict, optional
        Dictionary specifying column ranges to filter the DataFrame
        (e.g., {'col1': (min1, max1), 'col2': (min2, max2)})
    columns : list of int or str, optional
        Columns to select from each file
    prefix : str, optional
        Prefix to add to filenames when renaming columns

    Returns
    -------
    tuple of (list or None, list of pd.DataFrame)
        - List of selected column names/indices, or None if not determined
        - List of DataFrames containing the combined data from all files

    Raises
    ------
    ValueError
        If any file cannot be read or if selected columns differ between files
    TypeError
        If input types are invalid
    """
    logger = logging.getLogger("mylogger")

    # Type check and normalize input to list
    if isinstance(list_of_files, (str, dict)):
        list_of_files = [list_of_files]

    if not isinstance(list_of_files, list):
        raise TypeError(
            f"list_of_files must be str, dict, or list, got {type(list_of_files)}"
        )

    if len(list_of_files) == 0:
        raise ValueError("list_of_files cannot be empty")

    # Type check columns
    if columns is not None:
        if not isinstance(columns, list):
            raise TypeError(f"columns must be list or None, got {type(columns)}")
        if not all(isinstance(c, (int, str)) for c in columns):
            raise TypeError("columns must contain only int or str values")

    # Process files and validate column consistency
    list_of_df = []
    columns_selected = None

    for i, f in enumerate(list_of_files):
        # Validate file entry type
        if not isinstance(f, (str, dict)):
            raise TypeError(f"File entry {i} must be str or dict, got {type(f)}")

        if isinstance(f, dict) and "filename" not in f:
            raise ValueError(f"File entry {i} is a dict but missing 'filename' key")

        # Read and select columns
        wf = f.copy()
        if prefix:
            if isinstance(f, dict):
                wf["filename"] = prefix + "/" + f["filename"]
            else:
                wf = prefix + "/" + f

        # Read file into DataFrame
        df = read_file_to_df(wf, mask=mask)
        current_columns, df_selected = select_columns(columns, df)

        df_selected, columns_selected = process_dataframe_dimensions(
            df_selected, current_columns, logger
        )

        # Validate column consistency across files
        if columns_selected is None:
            columns_selected = current_columns
        elif columns_selected != current_columns:
            raise ValueError(
                f"Selected columns differ for file {i}. "
                f"Expected {columns_selected}, got {current_columns}"
            )

        list_of_df.append(df_selected)

    dataframes_selected = []

    # Combine dataframes with renamed columns
    if len(list_of_files) == 1:
        dataframes_selected = list_of_df[0]

    else:
        for dfs in zip(*list_of_df):
            renamed_dfs = []
            for df, f in zip(dfs, list_of_files):
                # Extract filename from dict or use string directly
                filename = f["filename"] if isinstance(f, dict) else f

                # Rename columns to filename
                df_renamed = df.copy()
                df_renamed.columns = [f"{c} @ {filename}" for c in df.columns]
                renamed_dfs.append(df_renamed)

            # Concatenate the renamed dataframes
            dataframes_selected.append(pd.concat(renamed_dfs, axis=1))

    if skim is not None:
        logger.debug(f"Skimmig data - keeping every {skim} value ...")
        for i in range(len(dataframes_selected)):
            dataframes_selected[i] = skim_dataframe(dataframes_selected[i], skim)

    return columns_selected, dataframes_selected


def read_file_to_df(
    fname: Union[str, Dict[str, str]], ftype: Optional[str] = None, **kwargs
) -> pd.DataFrame:
    """
    Read a file into a pandas DataFrame based on file type.

    Parameters
    ----------
    fname : str or dict
        Path to the file, or dict with 'filename' and 'type' keys
    ftype : str, optional
        File type/extension. If None, extracted from fname
    **kwargs : dict
        Additional arguments passed to pandas read functions

    Returns
    -------
    pd.DataFrame
        DataFrame containing the file data

    Raises
    ------
    ValueError
        If fname is None or file cannot be read
    """
    logger = logging.getLogger("mylogger")

    # Check filename
    if fname is None:
        logger.error("Filename must be provided")
        raise ValueError("Filename must be provided")

    # Handle dict input
    if isinstance(fname, dict):
        try:
            ftype = fname.get("type")
            fname = fname.get("filename")
        except Exception as e:
            logger.error(f"Invalid file descriptor dictionary: {e}")
            raise ValueError(f"Invalid file descriptor dictionary: {e}")

    # Log file reading
    logger.info(f"Reading data from {fname}...")

    # Extract file type from filename if not provided
    if ftype is None:
        parts = fname.split(".")
        ftype = parts[-1] if len(parts) > 1 else "out"
    else:
        ftype = ftype.lower()

    # If remote file, copy locally
    tmp_path = None
    if "@" in fname:
        # sp.run(["scp", fname, "."])
        # fname = fname.split("@")[-1].split(":")[-1].split("/")[-1]
        with tempfile.NamedTemporaryFile(delete=False) as tmp:
            tmp_path = tmp.name
            logger.debug(f"Temporary local file {tmp_path}")
        try:
            sp.run(["scp", fname, tmp_path])
            fname = tmp_path
        except Exception as e:
            raise ValueError(f"Failed to fetch remote file {fname}\n{str(e)}")

    # Check file existence
    if not os.path.exists(fname):
        msg = f"Input file {fname} does not exist"
        logger.error(msg)
        raise FileNotFoundError(msg)

    try:
        if ftype == "csv":
            remove_hash_from_first_line(fname)
            data = pd.read_csv(fname, sep=",", comment="#")

        elif ftype in ["lmp", "lammps"]:
            data = read_lammps_output_to_df(fname)

        elif ftype == "txt":
            data = pd.read_csv(
                fname,
                sep=r"\t+",
                comment="#",
                skip_blank_lines=True,
                header=None,
                engine="python",
            )
        else:
            # Default: whitespace-separated values
            data = pd.read_csv(
                fname,
                sep=r"\s+",
                comment="#",
                skip_blank_lines=True,
                header=None,
                engine="python",
            )

        # Handle multiple data blocks
        if isinstance(data, list):
            logger.warning(
                f"Multiple data blocks found in {fname}. Using the last block."
            )
            data = data[-1]

        # Log data info
        logger.info(f"... Number of rows: {data.shape[0]}")
        logger.info(f"... Number of columns: {data.shape[1]}")
        logger.debug(
            "Column selections must be 1-based integers or comma-separated integers."
        )
        msg = "Available columns:\n"
        msg += "\t0: Index\n"
        for i, col in enumerate(data.columns):
            msg += f"\t{i + 1}: {col}\n"
        logger.info(msg[:-1])  # Remove trailing newline

        # for l in data.columns:
        logger.debug(f"Data preview:\n{data.head()}")

        # Filter data by range if requested
        if kwargs.get("mask", None) is not None:
            data = mask_df_by_range(data, kwargs["mask"])

        # Dash formatting
        fmt = formatting()
        logger.info(fmt.dashes)

        return data

    except Exception as e:
        raise ValueError(f"Failed to read file {fname}: {str(e)}")

    finally:
        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)


def read_lammps_output_to_df(file_path: str) -> Union[pd.DataFrame, List[pd.DataFrame]]:
    """
    Read LAMMPS output file and parse data blocks into DataFrame(s).

    Parameters
    ----------
    file_path : str
        Path to the LAMMPS output file

    Returns
    -------
    pd.DataFrame or list of pd.DataFrame
        Single DataFrame if one block found, otherwise list of DataFrames

    Raises
    ------
    ValueError
        If no data blocks starting with 'Step' are found
    """
    logger = logging.getLogger("mylogger")
    logger.info("Reading LAMMPS output...")

    with open(file_path, "r") as file:
        lines = file.readlines()

    data_blocks = []
    i = 0

    while i < len(lines):
        if lines[i].strip().startswith("Step"):
            header_index = i

            # Find the end of this block
            footer_index = -1
            for j in range(header_index + 1, len(lines)):
                if lines[j].strip().startswith(("Loop", "WARNING")):
                    footer_index = j
                    break

            if footer_index == -1:
                data_lines = lines[header_index:]
                i = len(lines)
            else:
                data_lines = lines[header_index:footer_index]
                i = footer_index + 1

            # Convert extracted portion into DataFrame
            data_text = "".join(data_lines)
            try:
                df = pd.read_csv(io.StringIO(data_text), sep=r"\s+")
                data_blocks.append(df)
            except Exception as e:
                logger.warning(f"Failed to parse a data block: {e}")
        else:
            i += 1

    if not data_blocks:
        raise ValueError("Could not find any lines starting with 'Step'")

    logger.info(f"... Data blocks found: {len(data_blocks)}")

    return data_blocks[-1] if len(data_blocks) == 1 else data_blocks


def remove_hash_from_first_line(file_path: str) -> bool:
    """
    Remove '#' character from the first line of a file.

    Parameters
    ----------
    file_path : str
        Path to the file

    Returns
    -------
    bool
        True if successful, False otherwise
    """
    try:
        with open(file_path, "r") as file:
            lines = file.readlines()

        if not lines:
            return False

        lines[0] = lines[0].replace("#", "")

        with open(file_path, "w") as file:
            file.writelines(lines)

        return True
    except Exception as e:
        logging.getLogger("mylogger").error(f"Error removing hash: {e}")
        return False


def process_dataframe_dimensions(list_of_df, columns_selected, logger):
    """
    Process and reshape dataframes to ensure compatible dimensions.

    Args:
        list_of_df: List containing dataframe(s) or tuple of dataframes
        columns_selected: List of column selections to be updated
        logger: Logger instance for error reporting

    Returns:
        tuple: (renamed_dfs, updated_columns_selected)

    Raises:
        SystemExit: If more than 2 dimensions are provided
    """
    if not list_of_df or len(list_of_df) == 0:
        logger.error("Empty list_of_df provided")
        sys.exit(1)

    renamed_dfs = []
    # Extract dataframes based on structure
    if len(list_of_df) == 1:
        renamed_dfs.append(list_of_df[0].index.values)
        renamed_dfs.append(list_of_df[0])
    elif len(list_of_df) == 2:
        renamed_dfs.append(list_of_df[0])
        renamed_dfs.append(list_of_df[1])
    else:
        logger.error("More than 2 dimensions not supported for single file")
        sys.exit(1)

    # Validate dataframes
    if len(renamed_dfs) != 2:
        logger.error("Expected exactly 2 dataframes")
        sys.exit(1)

    for i, df in enumerate(renamed_dfs):
        if not hasattr(df, "shape"):
            logger.error(f"Element {i} is not a valid dataframe")
            sys.exit(1)

    # Handle dimension mismatch: replicate x_data if needed
    if renamed_dfs[0].shape[1] == 1 and renamed_dfs[1].shape[1] > 1:
        num_columns = renamed_dfs[1].shape[1]
        original_col_name = renamed_dfs[0].columns[0]

        # Replicate x_data to match y_data's columns
        renamed_dfs[0] = pd.concat([renamed_dfs[0]] * num_columns, axis=1)
        renamed_dfs[0].columns = [
            f"{original_col_name}_{i}" for i in range(num_columns)
        ]

        # Update columns_selected if it has entries
        if len(columns_selected) > 0 and len(columns_selected[0]) > 0:
            columns_selected[0] = [
                f"{columns_selected[0][0]}_{i}" for i in range(num_columns)
            ]

    # # Convert numeric column names to descriptive strings
    for i in range(len(columns_selected)):
        for j in range(len(columns_selected[i])):
            if is_numeric([columns_selected[i][j]]):
                columns_selected[i][j] = f"Column #{int(columns_selected[i][j])+1}"

    return renamed_dfs, columns_selected
