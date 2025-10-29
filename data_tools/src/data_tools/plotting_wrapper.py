from typing import Optional, Union, List
import numpy as np
import pandas as pd
import logging

from .plotting_with_plotly import _plot_plotly
from .plotting_with_matplotlib import _plot_matplotlib

# ============================================================================
# Generic Plotting Function
# ============================================================================


def plot(
    x: Union[
        list, np.ndarray, pd.DataFrame, List[Union[list, np.ndarray, pd.DataFrame]]
    ],
    y: Union[
        list, np.ndarray, pd.DataFrame, List[Union[list, np.ndarray, pd.DataFrame]]
    ],
    xlabel: str = "X",
    ylabel: str = "Y",
    title: Optional[str] = None,
    legend: Optional[Union[str, List[str]]] = None,
    logx: Optional[Union[bool, int]] = None,
    logy: Optional[Union[bool, int]] = None,
    grid: bool = False,
    ptype: str = "plotly",
    fout: Optional[str] = None,
    **kwargs,
) -> None:
    """
    Universal plotting function supporting both Plotly and Matplotlib.

    Parameters
    ----------
    x : array-like, DataFrame, or list of array-like/DataFrame
        X-axis data. Can be:
        - Single array: used for all y series
        - Single DataFrame: column(s) used for all y series
        - List of arrays/DataFrames: one per y series (must match y length)
    y : array-like, DataFrame, or list of array-like/DataFrame
        Y-axis data. Can be:
        - Single array: plots one series
        - Single DataFrame: plots one series per column
        - List of arrays/DataFrames: plots multiple series
    xlabel : str, default "X"
        X-axis label
    ylabel : str, default "Y"
        Y-axis label
    title : str, optional
        Plot title
    legend : str or list of str, optional
        Legend labels. If None and y contains DataFrames, uses column names.
        If list, must match number of y series
    logx : bool or int, optional
        Log scale for x-axis. If int, specifies base (default 10)
    logy : bool or int, optional
        Log scale for y-axis. If int, specifies base (default 10)
    grid : bool, default False
        Whether to show grid
    ptype : str, default "plotly"
        Plot type: "plotly" (interactive) or "matplotlib" (static)
    fout : str, optional
        Output filename. If None:
        - Plotly: shows in browser without saving
        - Matplotlib: saves to "plot.png"
    **kwargs : dict
        Additional styling arguments:
        - colors: color or list of colors
        - markers: marker style or list of marker styles
        - linestyles: line style or list of line styles (matplotlib)
        - customstyles: list of custom style codes ('l', 'lp', 'd', 'dp', 'p')
        - linewidth: line width or list of line widths (default 2)
        - markersize: marker size or list of marker sizes (default 15)

    Examples
    --------
    >>> # Single series
    >>> plot([1, 2, 3], [1, 4, 9])

    >>> # Multiple series with shared x
    >>> plot([1, 2, 3], [[1, 4, 9], [1, 2, 3]], legend=["Squared", "Linear"])

    >>> # DataFrame with automatic legend from column names
    >>> df = pd.DataFrame({'A': [1, 4, 9], 'B': [1, 2, 3]})
    >>> plot([1, 2, 3], df)  # Uses 'A' and 'B' as legend labels

    >>> # Multiple DataFrames
    >>> plot([df1, df2], [df3, df4])

    >>> # Custom line styles
    >>> plot([1, 2, 3], [[1, 4, 9], [1, 2, 3]], customstyles=['l', 'dp'])
    """
    logger = logging.getLogger("mylogger")
    logger.info(f"Creating {ptype} plot...")

    # Convert inputs to numpy arrays
    def to_array(data):
        if isinstance(data, pd.DataFrame):
            return data
        return np.array(data) if not isinstance(data, np.ndarray) else data

    def is_list_of_arrays(data):
        """Check if data is a list of arrays/lists/DataFrames"""
        if not isinstance(data, (list, tuple)) or len(data) == 0:
            return False
        first = data[0]
        if isinstance(first, str):
            return False
        if isinstance(first, (pd.DataFrame, pd.Series)):
            return True
        try:
            iter(first)
            return hasattr(first, "__len__")
        except (TypeError, AttributeError):
            return False

    def extract_series_from_data(data):
        """Extract individual series from data, handling DataFrames"""
        series_list = []
        column_names = []

        if is_list_of_arrays(data):
            for item in data:
                if isinstance(item, pd.DataFrame):
                    for col in item.columns:
                        series_list.append(item[col].values)
                        column_names.append(col)
                elif isinstance(item, pd.Series):
                    series_list.append(item.values)
                    column_names.append(item.name if item.name is not None else None)
                else:
                    series_list.append(to_array(item))
                    column_names.append(None)
        else:
            if isinstance(data, pd.DataFrame):
                for col in data.columns:
                    series_list.append(data[col].values)
                    column_names.append(col)
            elif isinstance(data, pd.Series):
                series_list.append(data.values)
                column_names.append(data.name if data.name is not None else None)
            else:
                series_list.append(to_array(data))
                column_names.append(None)

        return series_list, column_names

    # Extract y data and column names
    y_data, y_column_names = extract_series_from_data(y)
    n_series = len(y_data)

    # Handle x data
    if is_list_of_arrays(x):
        x_series, _ = extract_series_from_data(x)
        # If x has fewer series than y, repeat the last x for remaining y series
        if len(x_series) < n_series:
            if len(x_series) == 1:
                x_data = x_series * n_series
            else:
                raise ValueError(
                    f"Number of x series ({len(x_series)}) must be 1 or match number of y series ({n_series})"
                )
        elif len(x_series) > n_series:
            raise ValueError(
                f"Number of x series ({len(x_series)}) cannot exceed number of y series ({n_series})"
            )
        else:
            x_data = x_series
    else:
        x_series, _ = extract_series_from_data(x)
        if len(x_series) == 1:
            x_data = x_series * n_series
        else:
            # If x is a DataFrame with multiple columns, use each column for corresponding y
            if len(x_series) == n_series:
                x_data = x_series
            else:
                raise ValueError(
                    f"Number of x columns ({len(x_series)}) must be 1 or match number of y series ({n_series})"
                )

    # Handle legend labels
    if legend is None:
        labels = [None] * n_series
        # Use column names from y if available
        # for i, col_name in enumerate(y_column_names):
        #     if col_name is not None:
        #         labels.append(str(col_name))
        #     elif n_series > 1:
        #         labels.append(f"Series {i + 1}")
        #     else:
        #         labels.append(None)
    elif isinstance(legend, str):
        labels = [legend]
    else:
        labels = list(legend)
        if len(labels) != n_series:
            raise ValueError(
                f"Number of legend labels ({len(labels)}) must match number of series ({n_series})"
            )

    logger.info(f"Number of series to plot: {n_series}")

    # Validate data lengths
    logger.debug(f"x_data lengths: {[len(xi) for xi in x_data]}")
    logger.debug(f"y_data lengths: {[len(yi) for yi in y_data]}")
    for i in range(n_series):
        if len(x_data[i]) != len(y_data[i]):
            raise ValueError(
                f"Length of x and y data must match for series {i + 1} (got {len(x_data[i])} and {len(y_data[i])})"
            )
    logger.debug("Data validation passed.")

    # Route to appropriate plotting function
    if ptype.lower() == "plotly":
        _plot_plotly(
            x_data,
            y_data,
            labels,
            xlabel,
            ylabel,
            title,
            grid,
            logx,
            logy,
            fout,
            **kwargs,
        )

    elif ptype.lower() in ["matplotlib", "mpl"]:
        _plot_matplotlib(
            x_data,
            y_data,
            labels,
            xlabel,
            ylabel,
            title,
            grid,
            logx,
            logy,
            fout,
            **kwargs,
        )
    else:
        raise ValueError(f"Unknown plot type: {ptype}. Use 'plotly' or 'matplotlib'")
