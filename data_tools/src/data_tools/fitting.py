import numpy as np
import logging
import re
import sys
from lmfit.models import ExpressionModel
import itertools


def polynomial_expression(order):
    terms = [f"c{i}*x**{i}" for i in range(order + 1)]
    return " + ".join(terms), [f"c{i}" for i in range(order + 1)]


def find_good_init(model, x_data, y_data, param_names):
    best_params = None
    best_residual = np.inf

    values = [0.01, 0.1, 1.0, 20.0, 200.0, 1e4]
    values += [-0.01, -0.1, -1.0, -20.0, -200.0, -1e4]
    n = 3  # Change this to your desired length

    combinations = itertools.combinations_with_replacement(values, n)

    # Try multiple initializations
    for pset in combinations:
        params = model.make_params()
        for i, pname in enumerate(param_names):
            params[pname].set(value=pset[i], min=1e-10)

        try:
            result = model.fit(y_data, params, x=x_data)
            print(pset)
            print(best_residual)
            residual_sum = np.sum(result.residual**2)
            if residual_sum < best_residual and np.isfinite(residual_sum):
                best_residual = residual_sum
                best_params = params
        except:
            continue

    return best_params


def fitting_wrapper(X, Y, fit_type=None):

    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("mylogger")

    X = np.asarray(X)
    Y = np.asarray(Y)

    n_samples = Y.shape[0]

    if Y.ndim == 1:
        Y_columns = Y.reshape(-1, 1)
    else:
        Y_columns = Y.reshape(n_samples, -1)

    n_y_columns = Y_columns.shape[1]

    if X.ndim == 1:
        X_2d = X.reshape(-1, 1)
        X_columns = np.repeat(X_2d, n_y_columns, axis=1)
    else:
        X_2d = X.reshape(n_samples, -1)
        if X_2d.shape[1] == 1:
            X_columns = np.repeat(X_2d, n_y_columns, axis=1)
        elif X_2d.shape[1] == n_y_columns:
            X_columns = X_2d
        else:
            raise ValueError(
                f"X has {X_2d.shape[1]} columns but Y has {n_y_columns} columns. "
                "X must have 1 column or match Y's columns."
            )

    results = []

    poly_dict = {
        "linear": 1,
        "quadratic": 2,
        "parabola": 2,
        "cubic": 3,
        "quartic": 4,
        "quintic": 5,
    }

    custom_expressions = {
        "logarithmic": "a + b * log(x+c) + d",
        "logarithm": "a + b * log(x+c) + d",
        "log": "a + b * log(x+c) + d",
        "power_law": "a * (x+c)**b",
        "exponential": "a*exp(b*(x+c))",
        "exp": "a*exp(b*(x+c))",
        "arrhenius": "a*exp(-b/x)",
    }

    if isinstance(fit_type, int):
        expr_str, param_names = polynomial_expression(fit_type)

    elif isinstance(fit_type, str):
        fit_type = fit_type.strip().lower()

        if fit_type in list(poly_dict.keys()):
            order = poly_dict[fit_type]
            expr_str, param_names = polynomial_expression(order)

        elif fit_type in list(custom_expressions.keys()):
            expr_str = custom_expressions[fit_type]
            # Extract parameter names from expression
            param_names = sorted(
                set(re.findall(r"\b[a-zA-Z_]\w*\b", expr_str))
                - {"x", "np", "sin", "cos", "exp", "log", "sqrt"}
            )

        else:
            expr_str = fit_type
            param_names = sorted(
                set(re.findall(r"\b[a-zA-Z_]\w*\b", fit_type))
                - {"x", "np", "sin", "cos", "exp", "log", "sqrt"}
            )
    else:
        raise ValueError(
            "fit_type must be an int (polynomial degree) or a str (expression)."
        )

    expr_str = expr_str.replace("*x**0", "")
    expr_str = expr_str.replace("*x**1", "*x")

    for i in range(n_y_columns):
        x = X_columns[:, i]
        y = Y_columns[:, i]

        model = ExpressionModel(expr_str)
        for i, pname in enumerate(param_names):
            model.set_param_hint(pname, value=0.4, min=-np.inf, max=np.inf)
        params = model.make_params()

        if fit_type in ["logarithmic", "logarithm", "log"]:
            params["c"].set(min=-min(x) + 1e-5)

        if fit_type in ["power_law"]:
            params["b"].set(value=1.8)
            params["c"].set(min=-min(x) + 1e-5)

        if fit_type in ["arrhenius"]:
            params["a"].set(min=0)
            params["b"].set(min=0)

        if fit_type in ["exponential", "exp"]:
            params["a"].set(min=0)
            params["b"].set(min=0)
            params["b"].set(max=100 / max(np.abs(x)))
            params["c"].set(value=-(x.min() + x.max()) * 0.5)
            params["c"].set(max=-x.min() * 0.9)
            params["c"].set(min=-x.max() * 1.1)

        try:
            result = model.fit(y, x=x, params=params)
            fit_params = {
                k: {
                    "value": v.value,
                    "stderr": v.stderr,
                    "init_value": v.init_value,
                    "vary": v.vary,
                    "min": v.min,
                    "max": v.max,
                }
                for k, v in result.params.items()
            }

            # Replace parameter names with best-fit values (no uncertainty)
            evaluated_expr = expr_str
            for k, v in result.params.items():
                evaluated_expr = re.sub(rf"\b{k}\b", f"{v.value:.6g}", evaluated_expr)

            x_fit = x
            fitted_values = result.best_fit.tolist()
            # Compute fitted values on 100 evenly-spaced points
            # xmin, xmax = x.min() * 0.9, x.max() * 1.1
            # if fit_type == "logarithmic":
            #     xmin = max(min(x) - 1e-5, xmin)

            # x_fit = np.linspace(xmin, xmax, 100)
            # fitted_values = model.eval(result.params, x=x_fit).tolist()

            results.append(
                {
                    "type": (
                        "polynomial" if isinstance(fit_type, int) else "expression"
                    ),
                    "original_expression": expr_str,
                    "evaluated_expression": evaluated_expr,
                    "parameters": fit_params,
                    "residuals": result.residual.tolist(),
                    "fitted_range": x_fit,
                    "fitted_values": fitted_values,
                    "fit_report": result.fit_report(),
                }
            )
        except Exception as e:
            logger.error(f"Fitting failed for column {i}: {e}")
            results.append(
                {"type": "expression", "expression": expr_str, "error": str(e)}
            )
            sys.exit(1)

        # After your fit...
        logger.result("=" * 70)
        logger.result(f"Fit results for column {i}")
        logger.result("=" * 70)
        logger.result(f"Model Expression: {expr_str}")
        logger.result(f"Fitted Expression: {evaluated_expr}")
        logger.result("-" * 70)
        logger.result(
            f"{'Parameter':<12} {'Value':>12} {'Std Error':>12} {'Initial':>12}"
        )
        logger.result("-" * 70)

        for pname, param in result.params.items():
            stderr_str = f"{param.stderr:.6g}" if param.stderr is not None else "N/A"
            logger.result(
                f"{pname:<12} {param.value:>12.6g} {stderr_str:>12} {param.init_value:>12.6g}"
            )

        logger.result("-" * 70)
        logger.result(f"Chi-square: {result.chisqr:.6g}")
        logger.result(f"Reduced Chi-square: {result.redchi:.6g}")
        logger.debug(f"Akaike Info Criterion: {result.aic:.6g}")
        logger.debug(f"Bayesian Info Criterion: {result.bic:.6g}")
        logger.debug(f"Number of function evaluations: {result.nfev}")
        logger.debug(f"Number of data points: {result.ndata}")
        logger.debug(f"Number of variables: {result.nvarys}")
        logger.result(f"Success: {result.success}")
        logger.result("=" * 70)

    return results
