import numpy as np
import matplotlib.pyplot as plt

# Define the function f(x) = 1 / (1 + x^2)
def f(x):
    return 1 / (1 + x**2)

# Given x-values and divided differences for different n values
x_values_n2 = np.array([-5, 0, 5])
div_diff_n2 = np.array([0.0384615, 0.192308, -0.0384615])

x_values_n4 = np.array([-5, -2.5, 0, 2.5, 5])
div_diff_n4 = np.array([0.0384615, 0.0397878, 0.061008, -0.0265252, 0.00530504])

x_values_n6 = np.array([-5, -3.33333, -1.66667, 0, 1.66667, 3.33333, 5])
div_diff_n6 = np.array([0.0384615, 0.0264644, 0.0248454, 0.0149446, -0.0131699, 0.00420316, -0.000840633])

x_values_n8 = np.array([-5, -3.75, -2.5, -1.25, 0, 1.25, 2.5, 3.75, 5])
div_diff_n8 = np.array([0.0384615, 0.0223428, 0.013956, 0.0117043, 0.000674338, -0.00489646, 0.00243964, -0.000687223, 0.000137445])

# Define function to compute Newton's polynomial using divided differences
def newton_polynomial(x_vals, div_diffs, x):
    n = len(x_vals)
    result = div_diffs[0]
    term = 1.0
    for i in range(1, n):
        term *= (x - x_vals[i - 1])
        result += div_diffs[i] * term
    return result

# Define x points for plotting
x_plot = np.linspace(-5, 5, 1000)

# Compute polynomial interpolations for each n
poly_n2 = [newton_polynomial(x_values_n2, div_diff_n2, x) for x in x_plot]
poly_n4 = [newton_polynomial(x_values_n4, div_diff_n4, x) for x in x_plot]
poly_n6 = [newton_polynomial(x_values_n6, div_diff_n6, x) for x in x_plot]
poly_n8 = [newton_polynomial(x_values_n8, div_diff_n8, x) for x in x_plot]

# Plot the exact function and the polynomials
plt.figure(figsize=(10, 7))

# Plot exact function
plt.plot(x_plot, f(x_plot), 'k-', label='Exact function f(x)', linewidth=2)

# Plot interpolations
plt.plot(x_plot, poly_n2, 'r--', label='Polynomial n=2', linewidth=1)
plt.plot(x_plot, poly_n4, 'g--', label='Polynomial n=4', linewidth=1)
plt.plot(x_plot, poly_n6, 'b--', label='Polynomial n=6', linewidth=1)
plt.plot(x_plot, poly_n8, 'm--', label='Polynomial n=8', linewidth=1)

# Customize plot
plt.title('Runge Phenomenon for Polynomial Interpolation of f(x) = 1 / (1 + x^2)')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.grid(True)

# Show plot
plt.show()
