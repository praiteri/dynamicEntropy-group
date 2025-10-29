import sys
import logging
import numpy as np

import logging
from .logger import setup_logger
from .basic_parser import basic_parser
from .read_files import read_file_to_df
from .select_columns import select_columns
from .plot import plot


def integrate(xvals, yvals, xlower, xupper):
    """
    Integrate y over x from xlower to xupper.
    Use trapz to integrate over points closest to xlower, xupper.
    the +1 to idx_max is for numpy half-open indexing.
    """
    import numpy as np

    idx_min = np.argmin(np.abs(xvals - xlower))
    idx_max = np.argmin(np.abs(xvals - xupper)) + 1
    result = np.trapz(yvals[idx_min:idx_max], x=xvals[idx_min:idx_max])
    return result


def analytic_pairing_free_energy(x, pi4e0, charge, eps, kt, c=0):
    e = pi4e0 * charge / x / eps - kt * np.log(4.0 * np.pi * x**2) + c
    return e


def get_units(units):
    if units == "omm":
        u = {
            "energy": "kJ/mol",
            "distance": "nm",
            "c0": 0.6022,  # 1/(nm**3)
            "kB": 0.00831446261815324,  # kJ/(K mol)
            "pi4e0": 138.93545302827604,  # nm kJ/(mol e**2)
        }
    elif units == "lmp":
        u = {
            "energy": "eV",
            "distance": "A",
            "c0": 0.0006022,  # 1/(A**3)
            "kB": 0.000086173,  # eV / K
            "pi4e0": 14.399645,  # A eV/(e**2)
        }
    else:
        raise Exception("Units must be specified with +u omm/lmp")

    return u


def dielectric_constant(temp=300, water=None):
    if water is None:
        return ["exp", "spcfw", "amoeba"]
    elif water == "exp":
        eps = 249.355 - 0.7877 * temp + 0.0007192 * temp**2
    elif water == "spcfw":
        eps = 181.464 - 0.338253 * temp
    elif water == "amoeba":
        eps = 412.207 - 1.5853 * temp + 0.0017023 * temp**2
    else:
        raise Exception("Cannot compute dielectric constant for water model" + water)
    return eps


def association_constant(data, params):
    logger = logging.getLogger("mylogger")
    logger.info("Computing association constant from FES ..")

    if not isinstance(data, np.ndarray):
        data = data.to_numpy()

    #
    if isinstance(params["align"], str):
        try:
            parts = params["align"].split(":")
            assert len(parts) == 2
            params["align"] = [float(parts[0]), float(parts[1])]
        except Exception as e:
            logger.error(f"Invalid alignment distance specification: {e}")
            sys.exit(1)

    # Define input units
    # logger.debug(my.fmt["debug"].format("Units", params["units"]))
    if params["units"] == "omm":
        if params["align"] is None:
            params["align"] = [0.8, 1.0]
    elif params["units"] == "lmp":
        if params["align"] is None:
            params["align"] = [8.0, 10.0]
    units = get_units(params["units"])

    # Thermal energy
    kT = units["kB"] * params["temp"]

    # Dielectric constant
    if params["water"] in dielectric_constant():
        t = float(params["temp"])
        params["eps"] = dielectric_constant(t, params["water"])

    # Bjerrum length
    if params["rb"] is None:
        a = units["pi4e0"] * params["charge"] / params["eps"]
        params["rb"] = -0.5 * a / (units["kB"] * params["temp"])
        assert params["rb"] > 0.0, "Intergal limit must be specified for qq >= 0"

    logger.debug(f".. units : {params['units']}")
    logger.debug(f".. c0 ({units['distance']}) : {units['c0']}")
    logger.debug(f".. pi4e0 ({units['energy']}) : {units['pi4e0']}")
    logger.info(f".. charge product (e^2) : {params['charge']}")

    logger.info(f".. temperature (K) : {params['temp']}")
    logger.debug(f".. kB ({units['energy']}/K) : {units['kB']}")
    logger.debug(f".. kT ({units['energy']}) : {kT}")

    logger.info(f".. alignment distance ({units['distance']}) : {params['align']}")
    logger.info(f".. integral limit ({units['distance']}) : {params['rb']}")

    logger.info(f".. dielectic medium : {params['water']}")
    logger.info(f".... epsilon : {params['eps']}")

    analytic = analytic_pairing_free_energy(
        data[:, 0], units["pi4e0"], params["charge"], params["eps"], kT
    )

    mask = (data[:, 0] > params["align"][0]) & (data[:, 0] < params["align"][1])
    tail_d = data[mask]
    tail_a = analytic[mask]
    energyShift = np.mean(tail_a) - np.mean(tail_d[:, 1])
    data[:, 1] += energyShift
    logger.debug(f".. energy shift to align ({units['energy']}) : {energyShift}")

    e = np.exp(-data[:, 1] / kT)
    res = integrate(data[:, 0], e, data[0, 0], params["rb"]) * units["c0"]
    logger.result(f".. association constant ({units['energy']}) : {-kT * np.log(res)}")
    if params["units"] == "lmp":
        logger.result(f".. association constant (kJ/mol) : {-kT * np.log(res) * 96.485}")

    return analytic, energyShift


def main():
    """Main function for command-line association constant calculation"""
    parser = basic_parser(description="Association constant from FES")

    params = {
        "units": "omm",
        "charge": 0,
        "temp": 300,
        "water": "spcfw",
        "eps": 1,
        "align": None,
        "rb": None,
    }

    parser.add_argument(
        "-u",
        "--units",
        type=str,
        default=params["units"],
        choices=["omm", "lmp"],
        help="Units of input data: omm (kJ, nm) or lmp (eV, A). Default: omm",
    )
    parser.add_argument(
        "--charge",
        "--qq",
        type=float,
        default=params["charge"],
        help="Charge product (e^2) for the two interacting species. Default: 0",
    )
    parser.add_argument(
        "-t",
        "--temp",
        type=float,
        default=params["temp"],
        help="Temperature (K). Default: 300",
    )
    parser.add_argument(
        "--water",
        type=str,
        default=params["water"],
        choices=["exp", "spcfw", "amoeba", "Custom"],
        help="Water model for dielectric constant. Default: spcfw",
    )
    parser.add_argument(
        "--eps",
        type=float,
        default=params["eps"],
        help="Dielectric constant of the medium. If not specified, computed from water model.",
    )
    parser.add_argument(
        "--align",
        "--ra",
        type=str,
        default=params["align"],
        help="Alignment distance (x1:x2) in input units to align FES to analytical form. Default: 0.8:1.0 (omm) or 8.0:10.0 (lmp)",
    )
    parser.add_argument(
        "--rb",
        type=float,
        default=params["rb"],
        help="Upper limit of integral (in input units). Default: computed from charge product and dielectric constant for qq < 0. Must be specified for qq >= 0.",
    )

    # Parse arguments
    args, unknown = parser.parse_known_args()

    # Set up logging
    logger = setup_logger(args=args)
    logger.info("Computing association constant from FES")

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
    data_for_ac = np.stack((x_data, y_data), axis=-1)

    # Process association constant parameters
    for attr in params.keys():
        if getattr(args, attr) is not None:
            params[attr] = getattr(args, attr)

    analytic, dg = association_constant(data_for_ac, params)

    if ptype is not None:
        if params["units"] == "omm":
            xlabel = "Distance (nm)"
            ylabel = "Free Energy (kJ/mol)"
        else:
            xlabel = "Distance (A)"
            ylabel = "Free Energy (eV)"

        plot(
            x_data * 2,
            [y_data + dg, analytic],
            xlabel=xlabel,
            ylabel=ylabel,
            title=f"Plot of {args.input['filename']}",
            legend=["Data", "Analytic"],
            ptype=ptype,
            fout=fout,
            grid=True,
            logx=False,
            logy=False,
            customstyles=["l", "d"],
        )


if __name__ == "__main__":
    main()
