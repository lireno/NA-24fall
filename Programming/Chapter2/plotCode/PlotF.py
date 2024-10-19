import matplotlib.pyplot as plt
import numpy as np
import os

def read_points_from_file(filename):
    points = {}
    current_m = None
    
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith("For m = "):
                current_m = int(line.split('=')[1].strip()[:-1])
                points[current_m] = []
            elif line:
                x, y = map(float, line.split())
                points[current_m].append((x, y))
    return points

def plot_points(points):
    for m, point_set in points.items():
        plt.figure(figsize=(10, 6))
        
        x_values = [point[0] for point in point_set]
        y_values = [point[1] for point in point_set]
        x_values.extend(reversed([-point[0] for point in point_set]))
        y_values.extend(reversed(y_values))
        plt.plot(x_values, y_values, label=f'm = {m*2}',color = 'b', linestyle='--', marker='o', markersize=1)

        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title(f'Bezier Interpolated Points for m = {m*2}')
        plt.legend()
        plt.grid(True)
        plt.xlim([-2, 2]) 
        plt.ylim([-1.5, 2]) 
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

def f_new(x, y):
    return x**2 + (3/2 * y - np.sqrt(np.abs(x)))**2 - 3

# Create a grid of x and y values
x_vals = np.linspace(-2, 2, 400)
y_vals = np.linspace(-1.5, 2, 400)
X_new, Y_new = np.meshgrid(x_vals, y_vals)

# Evaluate the new function on the grid
Z_new = f_new(X_new, Y_new)

# Plot the contour where the new function equals 0 (this represents the curve)
plt.contour(X_new, Y_new, Z_new, levels=[0], colors='g')
plt.title(r'$x^2 + \left(\frac{3}{2}y - \sqrt{|x|}\right)^2 = 3$')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')

# Display the plot
plt.show()

if __name__ == "__main__":
    # Assuming pointsInF.txt is in the same directory as this script
    filename = os.path.join(os.path.dirname(__file__), "pointsInF.txt")
    
    if os.path.exists(filename):
        points = read_points_from_file(filename)
        plot_points(points)
    else:
        print(f"Error: {filename} does not exist.")
