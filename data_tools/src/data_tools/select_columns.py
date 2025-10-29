import sys
import logging

from typing import List, Tuple, Optional
import pandas as pd
import logging
import sys


def select_columns(
    columns_to_select: Optional[List[str]], data: pd.DataFrame
) -> Tuple[List[List[int]], List[pd.DataFrame]]:
    """Select specified columns from the data based on user input.

    This function processes column specifications (either by name or 1-based index),
    validates them against the provided DataFrame, and returns the selected data
    organized by dimension.

    Parameters
    ----------
    columns_to_select : Optional[List[str]]
        List of comma-separated column specifications. Each element represents
        a dimension and can contain column names or 1-based indices.
        Example: ["1,2,3", "Temperature,Pressure"]
    data : pd.DataFrame
        Source DataFrame containing the data to select from.

    Returns
    -------
    Tuple[List[List[int]], List[pd.DataFrame]]
        A tuple containing:
        - List of lists with 1-based column indices for each dimension
        - List of DataFrames, each containing the selected columns for that dimension

    Raises
    ------
    ValueError
        If a specified column name is not found in the DataFrame.
    SystemExit
        If columns_to_select is None or empty, or if column indices are out of range.

    Notes
    -----
    - Column indices in the input are expected to be 1-based (user-friendly)
    - Column indices in the return value are 1-based for consistency
    - Internal DataFrame operations use 0-based indexing (pandas convention)
    """
    logger = logging.getLogger("mylogger")
    logger.info("Selecting columns...")

    # Validate input: ensure columns are specified
    if columns_to_select is None:
        logger.error("No columns specified")
        sys.exit(1)

    logger.debug(f"Raw columns input: {columns_to_select}")
    dimensions: int = len(columns_to_select)

    # Validate input: ensure at least one dimension is specified
    if dimensions == 0:
        logger.error("Incorrect columns specification in --columns")
        sys.exit(1)

    # Create a copy to avoid modifying the original list
    columns_selection: List[str] = [spec for spec in columns_to_select]

    # Initialize storage for parsed column indices and resulting DataFrames
    columns_selected: List[List[int]] = []
    dataframes_selected: List[pd.DataFrame] = []

    # Parse each dimension's column specification
    for i, column_spec in enumerate(columns_selection):
        dimension_columns: List[int] = []

        # Split comma-separated values and process each item
        for item in column_spec.split(","):
            item = item.strip()

            # Case 1: Item is a numeric string (1-based index)
            if item.isdigit():
                dimension_columns.append(int(item))
            # Case 2: Item is a column name
            else:
                if item in data.columns:
                    # Convert 0-based pandas index to 1-based for consistency
                    column_index: int = data.columns.get_loc(item) + 1
                    dimension_columns.append(column_index)
                else:
                    raise ValueError(f"Invalid column name: {item}")

        columns_selected.append(dimension_columns)

    # Flatten the nested list to validate all column indices
    all_column_indices: List[int] = [
        idx for sublist in columns_selected for idx in sublist
    ]

    # Validate: ensure no column index exceeds the DataFrame's column count
    # Note: indices are 1-based, so we compare (idx - 1) with data.shape[1]
    if any(idx - 1 >= data.shape[1] for idx in all_column_indices):
        logger.error("Column index out of range.")
        sys.exit(1)

    # Create DataFrames for each dimension using the selected columns
    for i, col_indices in enumerate(columns_selected):
        # Convert 1-based indices to 0-based for DataFrame iloc indexing
        zero_based_indices: List[int] = [idx - 1 for idx in col_indices]

        # Extract the specified columns from the original DataFrame
        df_selected: pd.DataFrame = data.iloc[:, zero_based_indices]
        dataframes_selected.append(df_selected)

        logger.debug(f"Selected column indices for dimension {i}: {col_indices}")

    # Log preview of selected data for debugging
    for i, df in enumerate(dataframes_selected):
        logger.debug(f"Selected data preview - dimension {i}:\n{df.head()}")

    # Change the names in columns_selected from indices to actual column names
    for idx in range(len(columns_selected)):
        for jdx in range(len(columns_selected[idx])):
            columns_selected[idx][
                jdx
            ] = f"{data.columns[columns_selected[idx][jdx] - 1]}"

    return columns_selected, dataframes_selected


# def select_columns(columns_to_select, data):
#     """Select specified x and y columns from the data.

#     Parameters
#     ----------
#     columns_to_select : list
#         Parsed command-line arguments containing 'x' and 'y' attributes.
#     data : pd.DataFrame
#         DataFrame containing the data.

#     Returns
#     -------
#     tuple
#         Tuple containing:
#         - List of selected column indices for each dimension
#         - List of DataFrames, each containing the selected columns for that dimension
#     """
#     logger = logging.getLogger("mylogger")
#     logger.info("Selecting columns...")

#     if columns_to_select is None:
#         logger.error("No columns specified")
#         sys.exit(1)

#     logger.debug(f"Raw columns input: {columns_to_select}")
#     dimensions = len(columns_to_select)

#     if dimensions == 0:
#         logger.error("Incorrect columns specification in --columns")
#         sys.exit(1)

#     columns_selection = [x for x in columns_to_select]

#     # Parse column indices or names
#     columns_selected = [""] * len(columns_to_select)
#     dataframes_selected = []

#     for i, x_columns in enumerate(columns_selection):
#         columns_selected[i] = []
#         for item in x_columns.split(","):
#             item = item.strip()
#             if item.isdigit():
#                 columns_selected[i].append(int(item))
#             else:
#                 if item in data.columns:
#                     columns_selected[i].append(data.columns.get_loc(item) + 1)
#                 else:
#                     raise ValueError(f"Invalid x column: {item}")

#     # Validate column indices
#     if any(idx - 1 >= data.shape[1] for idx in sum(columns_selected, [])):
#         logger.error("Column index out of range.")
#         sys.exit(1)

#     # Create DataFrames for each dimension
#     for i, col_indices in enumerate(columns_selected):
#         # Convert 1-based indices to 0-based for DataFrame indexing
#         zero_based_indices = [idx - 1 for idx in col_indices]
#         df_selected = data.iloc[:, zero_based_indices]
#         dataframes_selected.append(df_selected)
#         logger.debug(f"Selected columns indices for dimension {i}: {col_indices}")

#     for i, df in enumerate(dataframes_selected):
#         logger.debug(f"Selected data preview - dimension {i}:\n{df.head()}")

#     return columns_selected, dataframes_selected
