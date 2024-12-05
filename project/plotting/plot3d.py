import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_spline_3d(file_path, output_path):
    # Load data
    x, y, z = [], [], []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:  # Ensure the line contains 3 values
                x.append(float(parts[0]))
                y.append(float(parts[1]))
                z.append(float(parts[2]))

    # Plot 3D spline
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, label="Spline Curve")
    ax.set_title(f"Spline Curve: {os.path.basename(file_path)}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend()
    plt.grid(True)

    # Save plot
    plt.savefig(output_path)
    plt.close()
    print(f"Saved 3D plot to {output_path}")

def process_all_splines_3d(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for file_name in os.listdir(input_dir):
        if file_name.endswith(".txt"):
            file_path = os.path.join(input_dir, file_name)
            output_path = os.path.join(output_dir, file_name.replace(".txt", ".png"))
            plot_spline_3d(file_path, output_path)

if __name__ == "__main__":
    input_directory = "./plotting/data/"
    output_directory = "./plotting/figures/"
    process_all_splines_3d(input_directory, output_directory)
