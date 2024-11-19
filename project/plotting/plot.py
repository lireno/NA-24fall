import os
import matplotlib.pyplot as plt

def plot_spline(file_path, output_path):
    # Load data
    x, y = [], []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                x.append(float(parts[0]))
                y.append(float(parts[1]))

    # Plot
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, label="Spline Curve")
    plt.scatter(x, y, color="red", s=10, label="Points")
    plt.title(f"Spline Curve: {os.path.basename(file_path)}")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.legend()
    plt.grid(True)
    plt.savefig(output_path)
    plt.close()
    print(f"Saved plot to {output_path}")

def process_all_splines(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for file_name in os.listdir(input_dir):
        if file_name.endswith(".txt"):
            file_path = os.path.join(input_dir, file_name)
            output_path = os.path.join(output_dir, file_name.replace(".txt", ".png"))
            plot_spline(file_path, output_path)

if __name__ == "__main__":
    input_directory = "./plotting/data/"
    output_directory = "./plotting/figures/"
    process_all_splines(input_directory, output_directory)
