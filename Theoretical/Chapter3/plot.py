import numpy as np
import matplotlib.pyplot as plt

def B2_i(x, i):
    # Define the knot values based on t_i = i
    t_i = i
    t_i_minus_1 = t_i - 1
    t_i_plus_1 = t_i + 1
    t_i_plus_2 = t_i + 2

    if t_i_minus_1 <= x < t_i:
        return ((x - t_i_minus_1) ** 2) / ((t_i_plus_1 - t_i_minus_1) * (t_i - t_i_minus_1))
    elif t_i <= x < t_i_plus_1:
        return ((x - t_i_minus_1) * (t_i_plus_1 - x)) / ((t_i_plus_1 - t_i) * (t_i_plus_1 - t_i_minus_1)) + \
               ((t_i_plus_2 - x) * (x - t_i)) / ((t_i_plus_2 - t_i) * (t_i_plus_1 - t_i))
    elif t_i_plus_1 <= x < t_i_plus_2:
        return ((t_i_plus_2 - x) ** 2) / ((t_i_plus_2 - t_i) * (t_i_plus_2 - t_i_plus_1))
    else:
        return 0

# Generate x values
x_values = np.linspace(-2, 6, 1000)

# Plot B^2_i for i = 0, 1, 2, 3
plt.figure(figsize=(10, 6))
for i in range(4):
    y_values = [B2_i(x, i) for x in x_values]
    plt.plot(x_values, y_values, label=f'$B^2_{i}$')

plt.xlabel('x')
plt.ylabel('Basis Functions')
plt.title('B-spline Basis Functions $B^2_i$ for i = 0, 1, 2, 3')
plt.ylim(0, 1.1)  # Set y-limits to avoid showing negative values
plt.grid()
plt.legend()
plt.show()
