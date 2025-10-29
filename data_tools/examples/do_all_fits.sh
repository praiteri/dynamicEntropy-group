python generate_random_data.py
qp --i function_fitting_data.csv --cols 1 2 --fit linear
qp --i function_fitting_data.csv --cols 1 3 --fit quadratic
qp --i function_fitting_data.csv --cols 1 4 --fit cubic
qp --i function_fitting_data.csv --cols 1 5 --fit power_law
qp --i function_fitting_data.csv --cols 1 6 --fit exponential
qp --i function_fitting_data.csv --cols 1 7 --fit logarithm
qp --i function_fitting_data.csv --cols 1 8 --fit arrhenius
