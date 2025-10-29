from typing import Optional, Union, List, Dict, Any

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

import logging

from .plotting_utilities import parse_custom_linestyle, _normalize_style_param


def set_plot_style_plotly(
    figsize: tuple = (1200, 800),
    fonttype: str = "serif",
    fontsize: int = 24,
    title_size: Optional[int] = None,
    label_size: Optional[int] = None,
    tick_size: Optional[int] = None,
    legend_size: Optional[int] = None,
    line_width: float = 2,
    marker_size: float = 15,
    grid: bool = True,
    template: str = "plotly_white",
) -> Dict[str, Any]:
    """
    Create a custom Plotly layout template for publication-quality plots.

    Parameters
    ----------
    figsize : tuple, default (1200, 800)
        Figure size in pixels (width, height)
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
        Width of plot lines and axis lines
    marker_size : float, default 15
        Size of plot markers
    grid : bool, default True
        Whether to display grid lines
    template : str, default "plotly_white"
        Base Plotly template ('plotly', 'plotly_white', 'plotly_dark', etc.)

    Returns
    -------
    dict
        Layout dictionary for Plotly figures
    """
    # Set default sizes
    title_size = title_size or int(fontsize * 1.2)
    label_size = label_size or fontsize
    tick_size = tick_size or fontsize
    legend_size = legend_size or int(fontsize * 0.8)

    # Configure font family
    if fonttype == "serif":
        font_family = "Computer Modern, serif"
    else:
        font_family = "Arial, sans-serif"

    layout = {
        "template": template,
        "width": figsize[0],
        "height": figsize[1],
        "font": dict(family=font_family, size=fontsize, color="black"),
        "title_font": dict(size=title_size),
        "xaxis": dict(
            title_font=dict(size=label_size),
            tickfont=dict(size=tick_size),
            showline=True,
            linewidth=line_width,
            linecolor="black",
            mirror=True,
            ticks="outside",
            tickwidth=line_width,
            showgrid=grid,
            gridwidth=line_width / 2,
            gridcolor="lightgray",
        ),
        "yaxis": dict(
            title_font=dict(size=label_size),
            tickfont=dict(size=tick_size),
            showline=True,
            linewidth=line_width,
            linecolor="black",
            mirror=True,
            ticks="outside",
            tickwidth=line_width,
            showgrid=grid,
            gridwidth=line_width / 2,
            gridcolor="lightgray",
        ),
        # "legend": dict(
        #     font=dict(size=legend_size),
        # ),
        "plot_bgcolor": "white",
        "paper_bgcolor": "white",
        "margin": dict(l=80, r=40, t=80, b=80),
    }

    # Set default trace properties for line width and marker size
    # These would need to be applied when creating traces
    # layout["_default_line_width"] = line_width
    # layout["_default_marker_size"] = marker_size

    return layout


def _plot_plotly(
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
    """
    Internal function for creating plots using Plotly.

    Parameters
    ----------
    x_data : list of array-like
        List containing x-axis data for each series.
    y_data : list of array-like
        List containing y-axis data for each series.
    labels : list of str
        List of labels for each data series.
    xlabel : str
        Label for the x-axis.
    ylabel : str
        Label for the y-axis.
    title : str
        Title of the plot.
    grid : bool
        Whether to show grid lines on the plot.
    logx : bool
        Whether to use logarithmic scale for the x-axis.
    logy : bool
        Whether to use logarithmic scale for the y-axis.
    fout : str or None
        Output file path. If None, the plot is shown in the browser.
        If ends with '.html', saves as HTML. Otherwise, saves as image.
    **kwargs : dict, optional
        Additional keyword arguments for customizing the plot:
            fontsize : int, default=24
                Font size for labels and title.
            linewidth : int or list, default=2
                Line width for each series.
            markersize : int or list, default=15
                Marker size for each series.
            customstyles : list of str, optional
                List of custom style codes for each series.
            colors : list, optional
                List of colors for each series.
            markers : list, optional
                List of marker symbols for each series.

    Returns
    -------
    None
        Displays the plot in the browser or saves it to a file.
    """
    logger = logging.getLogger("mylogger")
    n_series = len(y_data)

    # Extract kwargs with defaults
    fontsize = kwargs.get("fontsize", 24)
    linewidth = kwargs.get("linewidth", 2)
    markersize = kwargs.get("markersize", 15)
    customstyles = kwargs.get("customstyles", None)
    arragement = kwargs.get("arrangement", None)

    # Apply plotly template
    custom_layout = set_plot_style_plotly(
        fontsize=fontsize, line_width=linewidth, marker_size=markersize
    )

    # Default colors and markers
    default_colors = [
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
    ]
    default_markers = [
        "circle",
        "square",
        "diamond",
        "cross",
        "x",
        "triangle-up",
        "triangle-down",
        "star",
        "hexagram",
    ]

    # Handle customstyles if provided
    if customstyles:
        parsed_styles = [
            parse_custom_linestyle(code, "plotly") for code in customstyles
        ]
        logger.debug(f"Parsed custom styles: {parsed_styles}")

        dash_styles = [s["dash"] for s in parsed_styles]
        modes = [s["mode"] for s in parsed_styles]

        # Extract markers from modes - only use markers if mode includes 'markers'
        # For plotly, we derive marker presence from the mode, not a separate marker key
        style_markers = [
            (
                default_markers[i % len(default_markers)]
                if "markers" in s["mode"]
                else None
            )
            for i, s in enumerate(parsed_styles)
        ]

        # Use custom style markers if no explicit markers provided
        if "markers" not in kwargs:
            kwargs["markers"] = style_markers
    else:
        dash_styles = None
        modes = ["lines"] * len(
            y_data
        )  # Default mode; will be adjusted based on markers

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
    # For plotly: None means no marker symbol, but we need to check the mode
    # If mode includes 'markers', use default marker; otherwise keep as None
    if modes:
        markers = [
            (
                (m if m is not None else default_markers[i % len(default_markers)])
                if modes[i] and "markers" in modes[i]
                else None
            )
            for i, m in enumerate(markers)
        ]
    else:
        markers = [
            m if m is not None else default_markers[i % len(default_markers)]
            for i, m in enumerate(markers)
        ]

    linewidths = _normalize_style_param(linewidth, n_series, "linewidth")
    markersizes = _normalize_style_param(markersize, n_series, "markersize")

    logger.debug(f"Colors: {colors}")
    logger.debug(f"Markers: {markers}")
    logger.debug(f"Line widths: {linewidths}")
    logger.debug(f"Marker sizes: {markersizes}")

    # Create figure
    # fig = go.Figure()

    if arragement is None:
        nrows = 1
        ncols = 1
        arragement = [[1, 1] for i in range(len(y_data))]
    fig = make_subplots(
        rows=nrows,
        cols=ncols,
    )

    for i, (x, y, label) in enumerate(zip(x_data, y_data, labels)):
        color = colors[i]
        marker = markers[i]
        lw = linewidths[i]
        ms = markersizes[i]

        # Determine mode and dash style
        if modes is not None:
            mode = modes[i]
        else:
            mode = "lines" if marker is None else "lines+markers"

        if dash_styles is not None:
            dash = dash_styles[i]
        else:
            dash = "solid"

        # Configure marker and line based on mode
        marker_dict = None
        line_dict = None

        if "markers" in mode and marker is not None:
            marker_dict = dict(
                symbol=marker,
                size=ms,
                color=color,
                line=dict(width=1, color="white"),
            )

        if "lines" in mode:
            line_dict = dict(width=lw, color=color)
            if dash:
                line_dict["dash"] = dash

        if mode == "box":
            bin_size = x[1] - x[0]
            bar_width = bin_size * 0.8
            fig.add_trace(
                go.Bar(
                    x=x,
                    y=y,
                    name=label,
                    marker=dict(
                        color=color,  # transparent fill
                        line=line_dict,  # border/outline
                    ),
                    width=bar_width,
                    showlegend=(label is not None),
                ),
                row=1,
                col=1,
            )
        else:
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=y,
                    mode=mode,
                    name=label,
                    marker=marker_dict,
                    line=line_dict,
                    showlegend=(label is not None),
                ),
                row=arragement[i][0],
                col=arragement[i][1],
            )

    # Apply log scales
    xaxis_type = "log" if logx else "linear"
    yaxis_type = "log" if logy else "linear"

    # Update layout
    custom_layout["xaxis"]["showgrid"] = grid
    custom_layout["yaxis"]["showgrid"] = grid
    custom_layout["xaxis"]["type"] = xaxis_type
    custom_layout["yaxis"]["type"] = yaxis_type

    fig.update_layout(
        **custom_layout,
        title=title,
        xaxis_title=xlabel,
        yaxis_title=ylabel,
        showlegend=(len(y_data) > 1 or labels[0] is not None),
        legend=dict(
            x=0.98,
            y=0.98,
            xanchor="right",
            yanchor="top",
            bgcolor="rgba(255, 255, 255, 0.8)",
            bordercolor="black",
            borderwidth=1,
        ),
    )

    # Save or show
    if fout:
        if fout.endswith(".html"):
            fig.write_html(fout)
        else:
            fig.write_image(fout, width=1200, height=800, scale=2)
        logger.info(f"Plot saved to {fout}")
    else:
        pio.renderers.default = "browser"
        fig.show()
