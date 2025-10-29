import pandas as pd
import sys
import logging


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


def mask_df_by_range(df, params):
    if not isinstance(params, list):
        raise TypeError("params must be a list")

    for x in params:
        df = mask_df_by_range_w(df, x.split("="))

    return df


def mask_df_by_range_w(df, params):
    """
    Select rows from a DataFrame based on a column's values within a specified range.

    Parameters:
    -----------
    df : pandas.DataFrame
        The DataFrame to filter
    params : list
        A list containing [column_index, range_string] where:
            - column_index: int or string, the column to filter on
            - range_string: string in one of these formats:
                - "min:max:stride" (e.g., "1:20:2")
                - "min:max" (stride defaults to 1) (e.g., "1:20")
                - "value" (single value) (e.g., "10")
                - ":max" (from beginning to max) (e.g., ":12")
                - "min:" (from min to end) (e.g., "1:")
                - "::stride" (all values with stride) (e.g., "::2")
                - ":" (all values)

    Returns:
    --------
    pandas.DataFrame
        Filtered DataFrame containing only rows where the column values fall within the specified range
    """
    logger = logging.getLogger("mylogger")

    column_index, range_string = params
    logger.info(
        f"Selecting data by column `{column_index}` with range = [{range_string}]"
    )

    # Get the column values

    if is_numeric(column_index):
        column_index = int(column_index)

    if isinstance(column_index, int):
        if column_index > len(df.columns):
            logger.error(f"Column index out of bouds: {column_index}")
            sys.exit(1)

        if column_index == 0:
            column = pd.Series(df.index, name="index")
        else:
            column = df.iloc[:, column_index - 1]  # Convert 1-based to 0-based index
    else:
        try:
            column = df[column_index]
        except KeyError:
            logger.error(f"Unknown column name: {column_index}")
            sys.exit(1)

    # Parse the range string
    parts = range_string.split(":")

    # Initialize default values
    min_val = None
    max_val = None
    stride = 1

    # Parse different range formats
    if len(parts) == 1:
        # Single value format: "10"
        if parts[0] != "":
            min_val = max_val = float(parts[0])
    elif len(parts) == 2:
        # Min-max format: "1:20" or "1:" or ":20"
        if parts[0] != "":
            min_val = float(parts[0])
        if parts[1] != "":
            max_val = float(parts[1])
    elif len(parts) == 3:
        # Min-max-stride format: "1:20:2" or "::2"
        if parts[0] != "":
            min_val = float(parts[0])
        if parts[1] != "":
            max_val = float(parts[1])
        if parts[2] != "":
            stride = float(parts[2])

    logger.debug(f"Range parts: {parts}")
    logger.debug(f"Parsed range: min={min_val}, max={max_val}, stride={stride}")

    # Apply filtering based on parsed range
    mask = pd.Series(True, index=df.index)

    if min_val is not None:
        mask &= column >= min_val

    if max_val is not None:
        mask &= column <= max_val

    result = df[mask]
    if result.empty or result.shape[0] == 0:
        logger.error(
            f"Filtered DataFrame is empty. No rows match the criteria: {params}"
        )
        sys.exit(1)

    # Apply stride if specified
    if stride != 1:
        # Get indices of filtered dataframe
        indices = result.index.tolist()
        # Apply stride to indices
        strided_indices = (
            indices[:: int(stride)]
            if isinstance(stride, float) and stride.is_integer()
            else indices[:: int(stride)]
        )
        # Return dataframe with strided indices
        return df.loc[strided_indices]

    # Log the shape and content of the filtered DataFrame
    logger.debug(f"Filtered DataFrame shape: {result.shape}")
    logger.debug(f"Filtered DataFrame:\n{result}")

    return result
