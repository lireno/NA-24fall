import numpy as np
import matplotlib.pyplot as plt


days = [0, 6, 10, 13, 17, 20, 28]
sp1_weights = [6.67, 17.3, 42.7, 37.3, 30.1, 29.3, 28.7]
sp2_weights = [6.67, 16.1, 18.9, 15.0, 10.6, 9.44, 8.89]

# Divided differences provided
divided_differences_sp1 = np.array([6.67, 1.77167, 0.457833, -0.124778, 0.013566, -0.000978085, 4.1477e-05])
divided_differences_sp2 = np.array([6.67, 1.57167, -0.0871667, -0.0152729, 0.00257908, -0.000204804, 8.6768e-06])

# Function to calculate interpolation using divided differences
def interpolate(x, x_values, divided_differences):
    n = len(x_values)
    result = divided_differences[0]
    for i in range(1, n):
        term = divided_differences[i]
        for j in range(i):
            term *= (x - x_values[j])
        result += term
    return result

# Generate a range of x values for plotting the interpolation
x_range = np.linspace(0, 43, 1000)

# Calculate interpolated values for both species using divided differences
interp_sp1 = [interpolate(x, days, divided_differences_sp1) for x in x_range]
interp_sp2 = [interpolate(x, days, divided_differences_sp2) for x in x_range]

# Plot the larvae's weight data and interpolation for both species
plt.figure(figsize=(10, 6))

# Plot for Species 1
plt.plot(days, sp1_weights, 'bo', label='Sample 1 Data', color='blue')
plt.plot(x_range, interp_sp1, label='Sample 1 Interpolation', linestyle='--', color='blue')

# Plot for Species 2
plt.plot(days, sp2_weights, 'go', label='Sample 2 Data', color='green')
plt.plot(x_range, interp_sp2, label='Sample 2 Interpolation', linestyle='--', color='green')

# Add titles and labels
plt.title("Larvae's Weight in 43 days with Interpolated Functions")
plt.xlabel("Days")
plt.ylabel("Weight (g)")
plt.legend()
plt.grid(True)
plt.show()
