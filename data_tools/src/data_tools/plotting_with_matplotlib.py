from typing import Optional

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

import logging
from .plotting_utilities import parse_custom_linestyle, _normalize_style_param

# ============================================================================
# Plotting Utilities
# ============================================================================


def set_plot_style_matplotlib(
    figsize: tuple = (12, 8),
    fonttype: str = "serif",
    fontsize: int = 24,
    title_size: Optional[int] = None,
    label_size: Optional[int] = None,
    tick_size: Optional[int] = None,
    legend_size: Optional[int] = None,
    line_width: float = 2,
    marker_size: float = 15,
    grid: bool = True,
) -> None:
    """
    Set matplotlib rcParams for publication-quality plots.

    Parameters
    ----------
    figsize : tuple, default (12, 8)
        Figure size in inches (width, height)
    fonttype : str, default "serif"
        Font type to use ("serif" or other)
    fontsize : int, default 24
        Base font size
    title_size : int, optional
        Title font size (defaults to fontsize * 1.2)
    label_size : int, optional
        Axis label font size (defaults to fontsize)
    tick_size : int, optional
        Tick label font size (defaults to fontsize)
    legend_size : int, optional
        Legend font size (defaults to fontsize * 0.8)
    line_width : float, default 2
        Width of plot lines
    marker_size : float, default 15
        Size of plot markers
    grid : bool, default True
        Whether to display grid lines
    dpi : int, default 300
        Figure resolution
    """
    # Set default sizes
    title_size = title_size or int(fontsize * 1.2)
    label_size = label_size or fontsize
    tick_size = tick_size or fontsize
    legend_size = legend_size or int(fontsize * 0.8)

    # Configure serif font
    if fonttype == "serif":
        mpl.rcParams["font.family"] = "serif"
        cmfont = font_manager.FontProperties(
            fname=mpl.get_data_path() + "/fonts/ttf/cmr10.ttf"
        )
        mpl.rcParams["font.serif"] = cmfont.get_name()
        mpl.rcParams["mathtext.fontset"] = "cm"
        mpl.rcParams["axes.unicode_minus"] = False

    # Update all rcParams
    plt.rcParams.update(
        {
            "figure.figsize": figsize,
            # "figure.dpi": dpi,
            "axes.formatter.use_mathtext": True,
            "font.size": fontsize,
            "axes.titlesize": title_size,
            "axes.labelsize": label_size,
            "xtick.labelsize": tick_size,
            "ytick.labelsize": tick_size,
            "legend.fontsize": legend_size,
            "axes.linewidth": line_width,
            "lines.linewidth": line_width,
            "lines.markersize": marker_size,
            "xtick.major.width": line_width,
            "ytick.major.width": line_width,
            "xtick.minor.width": line_width,
            "ytick.minor.width": line_width,
            "xtick.major.size": 6.0,
            "ytick.major.size": 6.0,
            "xtick.minor.size": 3.0,
            "ytick.minor.size": 3.0,
            "axes.grid": grid,
            "grid.linestyle": "--",
            "grid.linewidth": line_width / 2,
            "grid.alpha": 0.7,
        }
    )


def _plot_matplotlib(
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
):
    """Internal function for Matplotlib plotting."""
    logger = logging.getLogger("mylogger")

    n_series = len(y_data)

    # Extract kwargs with defaults
    fontsize = kwargs.get("fontsize", 24)
    linewidth = kwargs.get("linewidth", 2)
    markersize = kwargs.get("markersize", 8)
    customstyles = kwargs.get("customstyles", None)
    arragement = kwargs.get("arrangement", None)

    set_plot_style_matplotlib(
        fontsize=fontsize, line_width=linewidth, marker_size=markersize
    )

    # Default colors
    default_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    default_markers = [
        "o",  # circle
        "s",  # square
        "D",  # diamond
        "^",  # triangle-up
        "h",  # hexagon (closest to hexagram)
        "v",  # triangle-down
        "p",  # pentagon
        "+",  # cross
        "x",  # x
        "*",  # star
    ]

    # Handle customstyles if provided
    if customstyles:
        parsed_styles = [
            parse_custom_linestyle(code, "matplotlib") for code in customstyles
        ]
        logger.debug(f"Parsed custom styles: {parsed_styles}")

        linestyles_from_custom = [s["linestyle"] for s in parsed_styles]
        markers_from_custom = [s["marker"] for s in parsed_styles]

        # Use custom style attributes if not explicitly provided
        if "linestyles" not in kwargs:
            kwargs["linestyles"] = linestyles_from_custom
        if "markers" not in kwargs:
            kwargs["markers"] = markers_from_custom

    # Normalize style parameters
    colors = _normalize_style_param(kwargs.get("colors", None), n_series, "colors")
    # Use default colors where None
    colors = [
        (
            default_colors[c]
            if isinstance(c, int)
            else (c if c is not None else default_colors[i % len(default_colors)])
        )
        for i, c in enumerate(colors)
    ]

    markers = _normalize_style_param(kwargs.get("markers", None), n_series, "markers")
    # For matplotlib: if marker is "o" (from customstyle), replace with default marker from cycle
    # If None, keep as None (no marker). If already a specific marker, keep it.
    markers = [
        default_markers[i % len(default_markers)] if m == "yes" else m
        for i, m in enumerate(markers)
    ]

    linestyles = _normalize_style_param(
        kwargs.get("linestyles", "-"), n_series, "linestyles"
    )

    linewidths = _normalize_style_param(linewidth, n_series, "linewidth")
    markersizes = _normalize_style_param(markersize, n_series, "markersize")

    logger.debug(f"Colors: {colors}")
    logger.debug(f"Markers: {markers}")
    logger.debug(f"Linestyles: {linestyles}")
    logger.debug(f"Line widths: {linewidths}")
    logger.debug(f"Marker sizes: {markersizes}")

    if arragement is None:
        nrows = 1
        ncols = 1
        arragement = [[0, 0] for i in range(len(y_data))]

    # Create figure
    fig, ax = plt.subplots(
        nrows,
        ncols,
        figsize=(12, 8),
    )

    if isinstance(ax, plt.Axes):
        # Single subplot case: convert to 2D array
        ax = np.array([[ax]])
    elif ax.ndim == 1:
        # 1D array case (1 row or 1 column): reshape to 2D
        ax = ax.reshape(1, -1) if ax.shape[0] > 1 else ax.reshape(-1, 1)

    for i, (x, y, label) in enumerate(zip(x_data, y_data, labels)):
        logger.debug(f"Plotting series {i+1}/{n_series}: label={label}")
        color = colors[i]
        marker = markers[i]
        linestyle = linestyles[i]
        lw = linewidths[i]
        ms = markersizes[i]

        # If no marker (None or empty string), set markersize to 0
        plot_ms = 0 if not marker else ms
        plot_marker = marker if marker else None

        idx = arragement[i][0]
        jdx = arragement[i][1]
        if marker == "box":
            bin_size = x[1] - x[0]
            bar_width = bin_size * 0.8
            ax[idx, jdx].bar(
                x,
                y,
                label=label,
                color=color,
                linestyle=linestyle,
                linewidth=lw,
                width=bar_width,
            )
        else:
            ax[idx, jdx].plot(
                x,
                y,
                label=label,
                color=color,
                marker=plot_marker,
                linestyle=linestyle,
                linewidth=lw,
                markersize=plot_ms,
            )

        # Set labels and title
        ax[idx, jdx].set_xlabel(xlabel)
        ax[idx, jdx].set_ylabel(ylabel)
        if title:
            ax[idx, jdx].set_title(title)

        # Set log scales
        if logx:
            ax[idx, jdx].set_xscale("log")
        if logy:
            ax[idx, jdx].set_yscale("log")

        # Grid
        if grid:
            ax[idx, jdx].grid(True, alpha=0.3)

        # Legend
        if len(y_data) > 1 or labels[0] is not None:
            ax[idx, jdx].legend()

    plt.tight_layout()

    if fout is None:
        plt.show()
    else:
        plt.savefig(fout, dpi=300, bbox_inches="tight")
        logger.info(f"Plot saved to {fout}")

    plt.close()
