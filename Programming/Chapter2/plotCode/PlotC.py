import numpy as np
import matplotlib.pyplot as plt

# Define the function for this problem: f(x) = 1 / (1 + 25x^2)
def f_chebyshev(x):
    return 1 / (1 + 25 * x**2)

# Given Chebyshev x-values and divided differences for different n values
x_values_n5 = np.array([0.951057, 0.587785, 6.12323e-17, -0.587785, -0.951057])
div_diff_n5 = np.array([0.0423501, -0.169057, 1.42548, 2.61208, 2.7465])

x_values_n10 = np.array([0.987688, 0.891007, 0.707107, 0.45399, 0.156434, -0.156434, -0.45399, -0.707107, -0.891007, -0.987688])
div_diff_n10 = np.array([0.0393884, -0.0887389, 0.189679, -0.534305, 2.11681, 8.28743, 11.9543, 10.3568, 5.51277, 8.9925e-16])

x_values_n15 = np.array([0.994522, 0.951057, 0.866025, 0.743145, 0.587785, 0.406737, 0.207912, 6.12323e-17, -0.207912, -0.406737, -0.587785, -0.743145, -0.866025, -0.951057, -0.994522])
div_diff_n15 = np.array([0.0388699, -0.0800675, 0.134961, -0.231498, 0.449158, -1.04753, 2.33944, 14.0283, 6.93539, -47.5358, -140.628, -235.756, -302.208, -331.791, -333.619])

x_values_n20 = np.array([0.996917, 0.97237, 0.92388, 0.85264, 0.760406, 0.649448, 0.522499, 0.382683, 0.233445, 0.0784591, -0.0784591, -0.233445, -0.382683, -0.522499, -0.649448, -0.760406, -0.85264, -0.92388, -0.97237, -0.996917])
div_diff_n20 = np.array([0.0386906, -0.0773136, 0.120771, -0.17893, 0.271509, -0.441458, 0.785471, -1.44565, 1.11365, 23.7955, -24.128, -331.452, -965.9, -1734.36, -2309.92, -2462.38, -2166.92, -1552.44, -788.326, 3.59221e-12])

# Reusing the Newton polynomial function to compute Chebyshev interpolation
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
x_plot = np.linspace(-1, 1, 1000)

# Compute polynomial interpolations for each n
poly_n5 = [newton_polynomial(x_values_n5, div_diff_n5, x) for x in x_plot]
poly_n10 = [newton_polynomial(x_values_n10, div_diff_n10, x) for x in x_plot]
poly_n15 = [newton_polynomial(x_values_n15, div_diff_n15, x) for x in x_plot]
poly_n20 = [newton_polynomial(x_values_n20, div_diff_n20, x) for x in x_plot]

# Plot the exact function and the Chebyshev interpolations
plt.figure(figsize=(10, 7))

# Plot exact function
plt.plot(x_plot, f_chebyshev(x_plot), 'k-', label='Exact function f(x)', linewidth=2)

# Plot interpolations
plt.plot(x_plot, poly_n5, 'r--', label='Chebyshev Polynomial n=5', linewidth=1)
plt.plot(x_plot, poly_n10, 'g--', label='Chebyshev Polynomial n=10', linewidth=1)
plt.plot(x_plot, poly_n15, 'b--', label='Chebyshev Polynomial n=15', linewidth=1)
plt.plot(x_plot, poly_n20, 'm--', label='Chebyshev Polynomial n=20', linewidth=1)

# Customize plot
plt.title('Chebyshev Interpolation of f(x) = 1 / (1 + 25x^2)')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.grid(True)

# Show plot
plt.show()
