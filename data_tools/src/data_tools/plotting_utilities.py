from typing import Dict


def parse_custom_linestyle(code: str, library: str = "matplotlib") -> Dict[str, str]:
    """
    Convert custom linestyle code to matplotlib or plotly style components.

    Parameters
    ----------
    code : str
        Style code: 'l' (line), 'lp' (line+points), 'd' (dashed),
        'dp' (dashed+points), 'p' (points only)
    library : str, optional
        Either 'matplotlib' or 'plotly' (default: 'matplotlib')

    Returns
    -------
    dict
        Dictionary with style parameters appropriate for the specified library.
        - matplotlib: {'linestyle': str, 'marker': str or None}
        - plotly: {'dash': str or None, 'mode': str}
    """
    styles_config = {
        "l": {
            "matplotlib": {"linestyle": "-", "marker": None},
            "plotly": {"dash": "solid", "mode": "lines"},
        },
        "lp": {
            "matplotlib": {"linestyle": "-", "marker": "yes"},
            "plotly": {"dash": "solid", "mode": "lines+markers"},
        },
        "d": {
            "matplotlib": {"linestyle": "--", "marker": None},
            "plotly": {"dash": "dash", "mode": "lines"},
        },
        "dp": {
            "matplotlib": {"linestyle": "--", "marker": "yes"},
            "plotly": {"dash": "dash", "mode": "lines+markers"},
        },
        "p": {
            "matplotlib": {"linestyle": "", "marker": "yes"},
            "plotly": {"dash": None, "mode": "markers"},
        },
        "b": {
            "matplotlib": {"linestyle": "", "marker": "box"},
            "plotly": {"dash": None, "mode": "box"},
        },
    }

    default_styles = {
        "matplotlib": {"linestyle": "-", "marker": None},
        "plotly": {"dash": "solid", "mode": "lines"},
    }

    style = styles_config.get(code, {})
    return style.get(library, default_styles[library])


def _normalize_style_param(param, n_series, param_name):
    """
    Normalize a style parameter to a list of length n_series.

    Parameters
    ----------
    param : scalar or list
        Style parameter (color, marker, linestyle, etc.)
    n_series : int
        Number of series to plot
    param_name : str
        Name of parameter (for error messages)

    Returns
    -------
    list
        List of parameter values, one per series
    """
    if param is None:
        return [None] * n_series

    # Check if it's a list/tuple
    if isinstance(param, (list, tuple)):
        if len(param) != n_series:
            raise ValueError(
                f"Length of {param_name} ({len(param)}) must match number of series ({n_series})"
            )
        return list(param)

    # Scalar value - replicate for all series
    return [param] * n_series
