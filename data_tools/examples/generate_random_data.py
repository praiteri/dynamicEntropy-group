import numpy as np
import pandas as pd

def add_noise(y, noise_level):
    """Add Gaussian noise to data"""
    noise = np.random.normal(0, noise_level * np.abs(y).mean(), len(y))
    return y + noise

# Set random seed for reproducibility
np.random.seed(42)

# Generate X values (temperature in Kelvin)
n_points = 50
xmin = 200
ymax = 400
X = np.linspace(210, 400, n_points)

# Noise level (as percentage of signal)
noise_level = 0.01

# Y1: Linear (y = a + bx)
Y1 = 5 + 3 * X
Y1 = add_noise(Y1, noise_level)

# Y2: Parabola (y = a + bx + cx^2)
Y2 = Y1 + 0.2 * (X-xmin)**2
Y2 = add_noise(Y2, noise_level)

# Y3: Cubic (y = a + bx + cx^2 + dx^3)
Y3 = Y2 - 0.001 * (X-xmin)**3
Y3 = add_noise(Y3, noise_level)

# Y4: Power law (y = a * x^b)
Y4 = 0.5 * (X-xmin+1)**1.8
Y4 = add_noise(Y4, noise_level)

# Y5: Exponential (y = a * exp(b*x))
Y5 = 10 * np.exp(0.03 * (X-xmin))
Y5 = add_noise(Y5, noise_level)

# Y6: Logarithmic (y = a + b*ln(x))
Y6 = 5 + 200. * np.log(X-xmin+1)
Y6 = add_noise(Y6, noise_level)

# Y7: Arrhenius (y = A * exp(-Ea/(R*T)))
# Common form: k = A * exp(-Ea/RT)
# Using Ea = 50000 J/mol, R = 8.314 J/(mol·K), A = 1e10
R = 8.314
Ea = 50000
A = 1e10
Y7 = A * np.exp(-Ea / (R * X))
Y7 = add_noise(Y7, noise_level)

# Create DataFrame
data = {
    'X': X,
    'Y1_Linear': Y1,
    'Y2_Parabola': Y2,
    'Y3_Cubic': Y3,
    'Y4_Power': Y4,
    'Y5_Exponential': Y5,
    'Y6_Logarithmic': Y6,
    'Y7_Arrhenius': Y7
}

df = pd.DataFrame(data)

# Save to CSV
filename = 'function_fitting_data.csv'
df.to_csv(filename, index=False)

print(f"CSV file '{filename}' created successfully!")
print(f"\nData shape: {df.shape}")
print(f"\nFirst few rows:")
print(df.head())
print(f"\nColumn descriptions:")
print("X: Temperature values (200-400 K)")
print("Y1_Linear: Linear relationship")
print("Y2_Parabola: Quadratic relationship")
print("Y3_Cubic: Cubic polynomial")
print("Y4_Power: Power law relationship")
print("Y5_Exponential: Exponential growth")
print("Y6_Logarithmic: Logarithmic relationship")
print("Y7_Arrhenius: Arrhenius equation (reaction rate vs temperature)")
