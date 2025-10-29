import numpy as np
import logging
import sys

from lmfit import Model  # Curve fitting library
from scipy import integrate

from .logger import setup_logger
from .basic_parser import basic_parser
from .read_files import read_file_to_df
from .select_columns import select_columns
from .plot import plot


def tail(x, l=1e-5, c=-40):
    return (4 * np.pi * x**3 * l) / 3 + c


def kb_integral(data):
    """ """
    logger = logging.getLogger("mylogger")

    logger.info("Computing Kirkwood-Buff integral ..")

    if not isinstance(data, np.ndarray):
        data = data.to_numpy()

    first_nonzero = np.argmax(data[:, 1] > 0)
    x = data[first_nonzero:, 0]
    y = data[first_nonzero:, 1]

    # KB integral from the first non-zero value
    KB0 = -4 * np.pi * x[0] ** 3 / 3

    # limit model to fit the tail of the KB integral
    xlim = x[-1] * 0.6
    lim = np.argmax(x > xlim)

    # The tail is considered to start from 60% of the maximum value of x
    lmodel = Model(tail)
    params = lmodel.make_params()

    # Self consistent loop to find the shift to correct the g(r)
    shift = 0
    for i in range(100):
        intrgrand = 4 * np.pi * x**2 * (np.array(y) - 1 - shift)
        kb = 0.6022 * (integrate.cumulative_trapezoid(intrgrand, x) + KB0)
        result = lmodel.fit(kb[lim:], params, x=x[lim:-1])
        ds = result.params["l"].value
        if np.abs(ds) < 1e-8:
            break
        shift += ds

    logger.debug(f".. number of cycles: {i}")
    logger.debug(f".. final g(r) shift: {shift}")
    logger.result(f".. KB integral (cm^3/mol) : {result.params['c'].value}")

    return kb, result.params["c"].value


def main():
    """Main function for command-line association constant calculation"""
    parser = basic_parser(description="Association constant from FES")

    # Parse arguments
    args, unknown = parser.parse_known_args()

    # Set up logging
    logger = setup_logger(args=args)
    logger.info("Computing Kirkwood-Buff integral from g(r)")

    # Select plot type and file output for matplotlib
    fout = args.output if args.output else "plot.png"
    ptype = None
    if args.plottype in ["pl", "plotly"]:
        ptype = "plotly"
    elif args.plottype in ["pp", "pyplot", "matplotlib"]:
        ptype = "matplotlib"

    # Read data
    data = read_file_to_df(args.input, mask=args.select)

    # Parse column indices or names
    columns_selected, df_selected = select_columns(args.columns, data)

    # Validate column selections for plotting
    if len(columns_selected) != 2:
        logger.error("For association constant, provide one x and one y column.")
        sys.exit(1)

    # Prepare x and y data
    x_data = df_selected[0].iloc[:, 0].values
    y_data = df_selected[1].iloc[:, 0].values
    data_for_kb = np.stack((x_data, y_data), axis=-1)

    kb_func, kb_val = kb_integral(data_for_kb)

    n = len(kb_func)

    if ptype is not None:
        plot(
            [x_data[-n:]] * 2,
            [kb_func, [kb_val] * n],
            xlabel="r",
            ylabel=r"$\int 4\pi r^2 [g(r)-1]\mathrm{d}r$",
            title=f"Plot of {args.input['filename']}",
            legend=[None, f"KB = {kb_val:.3f} (cm^3/mol)"],
            ptype=ptype,
            fout=fout,
            grid=True,
            logx=False,
            logy=False,
            customstyles=["l", "d"],
        )


if __name__ == "__main__":
    main()
